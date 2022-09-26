% 6/28/2011: Changed calculation of bed slope. No longer consider beds
% flatter than 1/1000 or that get deeper close to shore to be flat. No
% longer 'smooth beds a bit.'

%% Set up %%%%%%%%%%%
clear all;clc;clf;
homedir = '/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model/analysis/Fit_Cd 110622';
%homedir = 'P:\Users\Malin Pinsky\NatCap\Fit_Cd 110622'
cdata = homedir;
global g rho Cf; g=9.81; rho=1024;

% Load data
infile='QuantforGreg_2011-07-04.txt';
cd(cdata);data=load(infile);cd(homedir);
len=length(data);tic

%Init. output vectrs: nHfit (number of wave measurements used), Cdval2(whole bed), Difval2(RMS err whole bed),R2val2 (R2 = 1-SSerror/SStotal),slo(bed slope),flatbed(check if transect flat)
% Bval2 (breaking coefficient)
Cdval2=-999*ones(len,1); slo=Cdval2;flatbed=Cdval2;
nHfit=Cdval2;Bval2=Cdval2;Difval2=Cdval2;R2val2=Cdval2;

%% Set up and fit models for each case

%for kk=1:len % step through each line of the input file
for kk=592:594

% a parallel loop version (on a multi-core computer)
%matlabpool open local 3
%parfor kk = 1:len

%% Initialize Vectors and default values
kk
% fig=1;if fig;clf;end

%% Bathymetry
xd=[0 data(kk,12:3:30)];xd(isnan(xd))=[];  %Location of Wave Height measurements (in m)
hd=data(kk,10:3:31);hd(isnan(hd))=[]; % depth at measurement locations

% Create x-axis (horiz scale depends on precision of distance given)
if(sum(mod(xd,0.5)) == 0); dx=0.5;  % if measurements fall on half-meters
elseif(sum(mod(xd,0.1)) == 0); dx=0.1; % if every 10 cm
elseif(sum(mod(xd,0.02)) == 0); dx=0.01; % if every 2 cm (e.g. Bouma et al. 2005)
elseif(sum(mod(xd,0.01)) == 0); dx=0.01; % use every 1 cm as a last resort
end
x=0:dx:xd(end); % Create an X-axis at the resolution we just determined
lx=length(x); L=length(xd);

%Check bed slope and interplate bed depths
p=polyfit(xd,hd,1);slo(kk)=1/(-p(1));  %Slope of bed and linear fit for bathy
if all(hd(1) == hd); flatbed(kk)=1; % label if bed is truly flat
else flatbed(kk)=0; %If not, bed is plane slope or something more complicated
end;
h=interp1(xd,hd,x);

%Plot bathy to check how meast and final profiles differ
% plot(x,-h,'linewidth',2);grid on;hold on;plot(xd,-hd,'or','linewidth',2);

%% Wave info
T=data(kk,4);sig=2*pi/T;fp=1/T;% Wave period, frequency etc. 
Hmeas=data(kk,8:3:29);Hmeas(isnan(Hmeas))=[]; % Wave height measured
if data(kk,3);irreg=1;else irreg=0;end %Indicates whether waves are (ir)regular. Matters to determine coeff. for transformation

%% Vegetation field
pos=xd*0;
for xx=1:length(xd)
    [~, pos(xx)]=closest(x,xd(xx));
end; %ID loc of H mst
[~,veg]=closest(x,data(kk,9)); %Loc of start of veg (as an index into x)

N=h*0; bv=h*0; hv=h*0;
x1=veg; x2=pos(end);%Boundaries of vegetation.  Assumes that all mst after index of 'veg' are taken in vegetated areas
N(x1:x2)=data(kk,5); %Shoot density in shoots/m2
bv(x1:x2)=data(kk,6); %Cross-sectional area of a single shoot in m2
hv(x1:x2)=data(kk,7); %Shoot height in m

%Remove first measurement that is outside of vegetation field, but only if
%there are TWO outside
out = find(xd < x(x1)); % find the index of wave measurements that are outside the vegetation
if length(out)>1;%If there are points that are before the veg starts
    out = out(1:(length(out)-1)); % remove the index closest to the veg
    newx1 = pos(max(out)+1); % index into x for new start:first H measurement inside veg
    N=N(newx1:x2);bv=bv(newx1:x2);hv=hv(newx1:x2);
    Hmeas(out)=[];xd(out)=[];hd(out)=[]; % wave heights (measured)
    xd=xd-xd(1); % location of wave height measurements (in meters)
    h=h(newx1:x2);x=x(newx1:x2)-x(newx1); % depth and distance vectors
    lx=length(x); L=length(xd);
end;
% or if the study is Bradley & Houser 2009 (kk=592:594), remove first
% measurement because of reflection effect
if ((kk == 592) || (kk==593) || (kk==594));
    out = 1;
    newx1 = pos(max(out)+1); % index into x for new start:first H measurement inside veg
    N=N(newx1:x2);bv=bv(newx1:x2);hv=hv(newx1:x2);
    Hmeas(out)=[];xd(out)=[];hd(out)=[]; % wave heights (measured)
    xd=xd-xd(1); % location of wave height measurements (in meters)
    h=h(newx1:x2);x=x(newx1:x2)-x(newx1); % depth and distance vectors
    lx=length(x); L=length(xd);    
end

h=smooth(h,round(lx/5)); if size(h,1)~=size(x,1);h=h';end; %Smooth h vector
pos=xd*0;for xx=1:L;[~,pos(xx)]=closest(x,xd(xx));end; %ID loc of H meast
            
%% Estimate Drag Coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(Hmeas)>1
    %% ID Location in depth vector where bed is flat and sloping
    %segs gives loc where sloping bed starts and ends (row vector)
    %segf gives loc where flat bed starts and ends (row vector)
    
    Dif=diff(hd); flat=find(~Dif);flat=union(flat,flat+1);xflat=[]; %ID where bed is flat is flagged by 1's
    if isempty(flat); flat=0;segf=0;segs=[1 L];Lf=0; %Bed has no flat portions
    else
        temp=flat*0;Lf=length(flat);
        for ii=1:Lf-1
            if flat(ii+1)~=flat(ii)+1;temp(ii)=1;end
        end; segf=find(temp);  %Locate where brk from 1 to 0 exists in the flat vector
        if isempty(segf);
            segf=[1 Lf];
        else a=[1 segf+1];b=[segf Lf];
            segf=[a' b']; segf=flat(segf); %segf gives loc where flat bed starts and ends (row vector)
        end;
        if mean(segf);xflat=xd(segf);end % get actual x-values for the beginning and end of the flat portion
        
        temp=setdiff(1:L,flat);
        slopy=union(temp,temp+1); slopy=union(slopy,temp-1);
        slopy(slopy<1)=[];slopy(slopy>L)=[]; %slopy gives loc of sloping bed
        if ~isempty(slopy);temp=slopy*0;Ls=length(slopy);
            for ii=1:Ls-1
                if slopy(ii+1)~=slopy(ii)+1;temp(ii)=1;end
            end; segs=find(temp);  %Locate where brk from 1 to 0 exists
            a=[1 segs+1];b=[segs Ls];
            segs=[a' b']; segs=slopy(segs);%segs gives loc where sloping bed starts and ends (row vector)
        else segs=0;
        end;
    end
    
    %% Transform waves over whole bed
    if size(x,1)~=1;x=x';end;
    temp=[];up=[];

    %Determine location of flat portion
    if segf~=0;
        for ss=1:size(segf,1); temp=[temp xd(segf(ss,1)):dx:xd(segf(ss,2))]; end
        xflat2=temp(:);
        % Why did Greg have these next two lines in? It causes any data
        % with dx = 0.01 to fail (since it rounds to nearest 0.1)
        %xflat2=round(xflat2*10)/10;
        %x=round(x*10)/10;
    else xflat2=1;
    end;%tic
        
    % Find Cd that minimize RMSE for each of B = 0.4 to 1.6
    Cf=.003; % Bottom friction coefficient
    Bmat=.4:0.1:1.6;%Breaking dissipation coefficients to try
    Cd2=Bmat*0;Dif2=Cd2; % initialize vectors to hold Cd and RMSE values for each value of B
    for bb=1:length(Bmat);%Iterate for each value of B
        B=Bmat(bb)
        f=@(Cd) Wave_Baldock_GamKnown2(Hmeas,fp,Cd,bv,N,hv,B,h,hd,flat,x,xflat,xflat2,pos);
        [Cdtemp1,Diftemp1] = fminsearch(f,100,optimset('TolX',1e-10, 'TolFun', 1e-10));
        [Cdtemp2,Diftemp2] = fminsearch(f,1,optimset('TolX',1e-10, 'TolFun', 1e-10));
        [Cdtemp3,Diftemp3] = fminsearch(f,0,optimset('TolX',1e-10, 'TolFun', 1e-10));
         
        Diftemp = [Diftemp1 Diftemp2 Diftemp3];
        inds = find(Diftemp == min(Diftemp)); % get all starting points that produced the minimum RMSE
        Cdtemp = [Cdtemp1 Cdtemp2 Cdtemp3];
        if(isempty(inds)) % set to NA if no minimum (e.g., all Diftemp1-3 are NaN)
            Cd2(bb) = -999;
            Dif2(bb) = NaN;
        else
            % make sure all min RMSEs suggest nearly the same Cd
            if(all(Cdtemp(inds) - mean(Cdtemp(inds)) < 1e-4) || all(Cdtemp(inds) > 1000))
                Cd2(bb) = mean(Cdtemp(inds));
                Dif2(bb) = mean(Diftemp(inds));
            else % otherwise, set to default NA value (didn't converge)
                Cd2(bb) = -999;
                Dif2(bb) = NaN;
            end
            
        end
    end
    
    % Find value of B where Cd is > 0 and that minimizes RMSE
    Dif2(Cd2<0)=999; % exclude negative Cd values
    [m1,r]=min(Dif2); % find the minimum RMSE
    Bval2(kk)=Bmat(r); % save the best B value from this line of data
    Cdval2(kk)=Cd2(r); % save the best Cd value
    Difval2(kk)=Dif2(r); % save the best RMSE value
    R2val2(kk)=1-((Dif2(r)^2)*L)/sum((Hmeas-mean(Hmeas)).^2); % R2 = 1-SSerror/SStotal
    
    % Plot observed and fit
    % hcalc =
    % Wave_Baldock_GamKnown_H2(Hmeas(1),fp,Cdval2(kk),bv,N,hv,Bval2(kk),h,hd,flat,x,xflat,xflat2);
    % plot(x, -h, 'blue'); hold on; 
    % plot(xd, Hmeas, '*'); hold on;
    % plot(x, hcalc,'red')

else Cdval2(kk)=-999; Bval2(kk)=-999;Difval2(kk)=-999;R2val2=-999;% if one or fewer measured wave heights
end

nHfit(kk) = L; % number of wave measurements used for this fit 

end % final end for the loop through each line of data file
toc;

% if using a multicore version
%matlabpool close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute KC and Re #
if 1
    % Load data
    %cd(cdata);data = load('VegData2.txt');cd(homedir); % not needed, since already loaded earlier
    len=length(data);Re=zeros(8,len)+NaN;KC=Re;Re2=Re;KC2=Re; % Re is Reynolds # from stem diameter, KC is KC #, Re2 is Reynodsl # from stem diamter x density
    ubm = zeros(1,len); ksave = ubm; sigsave = ubm; Hsave = ubm;hdsave = ubm;
    
    for kk=1:len
        hd=data(kk,10:3:31);hd(isnan(hd))=[]; hdsave(kk) = hd(1); % depth at measurement locations
        A= data(kk,6); %Diameter in m
        D = data(kk,5); % Density in stems/m2
        hv=data(kk,7); % Height of vegetation
        
        
        T = data(kk,4); sig=2*pi/T; % Wave period
        Hmeas=data(kk,8:3:29);Hmeas(isnan(Hmeas))=[];Hsave(kk)=Hmeas(1);
        x1=data(kk,9);%indicates if meast outside veg
        if x1>0;Hmeas(1)=[];hd(1)=[];end % trim off first measurement if it's outside the vegetation (won't trim off second even if it's also outside veg)
        
        for xx=1:length(hd)
            k=iterativek(sig,hd(xx)); ksave(kk) = k; sigsave(kk)=sig;
            hdhv = hd(xx)-hv; hdhv = max([hdhv 0]); % depth can't be < 0
            ubm(kk)=0.5*Hmeas(xx)*sig*cosh(k*hdhv)/sinh(k.*hd(xx));%Wave orb. velo. at top of veg
            Re(xx,kk)=A.*ubm(kk)/1.17e-6;
            Re2(xx,kk)=(A*D).*ubm(kk)/1.17e-6;

            ubmKC=0.5*Hmeas(xx)*sig*1/sinh(k.*hd(xx));%Max wave orb. velo.
            KC(xx,kk)=T*ubmKC./A;
            KC2(xx,kk)=T*ubm(xx)./A;
        end
    end
end
% figure(2);
% subplot 311;semilogy(1:len,Re,'o',1:len,nanmean(Re),'b','linewidth',2);grid on;
% ylabel('Reynolds #','fontsize',12,'fontweight','bold');
% subplot 312;semilogy(1:len,KC,'d',1:len,nanmean(KC),'b','linewidth',2);grid on;
% ylabel('K-C #','fontsize',12,'fontweight','bold');
% xlabel('Case #','fontsize',12,'fontweight','bold')
% subplot 313;semilogy(1:len,Re2,'d',1:len,nanmean(Re2),'b','linewidth',2);grid on;
% ylabel('Re2 #','fontsize',12,'fontweight','bold');
% xlabel('Case #','fontsize',12,'fontweight','bold')


ReMean = -999*ones(len,1); ReSd = ReMean;KCMean = ReMean; KCSd = ReMean;Re2Mean = ReMean; Re2Sd = ReMean;KC2Mean = ReMean; KC2Sd = ReMean;
for kk = 1:len
    temp = Re(:,kk); temp(isnan(temp))=[];
    ReMean(kk) = mean(temp);
    ReSd(kk) = std(temp);

    temp = KC(:,kk); temp(isnan(temp))=[];
    KCMean(kk) = mean(temp);
    KCSd(kk) = std(temp);

    temp = Re2(:,kk); temp(isnan(temp))=[];
    Re2Mean(kk) = mean(temp);
    Re2Sd(kk) = std(temp);

    temp = KC2(:,kk); temp(isnan(temp))=[];
    KC2Mean(kk) = mean(temp);
    KC2Sd(kk) = std(temp);
end



%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
% append new data columns to the data that were read in
out = [data, nHfit, Cdval2, Difval2, R2val2, Bval2, slo, flatbed, ReMean, ReSd, KCMean, KCSd, KC2Mean, KC2Sd]; 
cd(cdata);

% make a header of column names
header = fgetl(fopen(infile)); fclose('all');
header = regexprep(header, '\t', ',');
header = [header, ',nHfit,Cdval2,Difval2,R2val2,Bval2,slo,flatbed,ReMean,ReSd,KCMean,KCSd,KC2Mean,KC2Sd'];

% write out a csv
dlmwrite(['VegData2Out_', date(), '.csv'], header, 'delimiter', '', 'precision', 10)
dlmwrite(['VegData2Out_', date(), '.csv'], out, '-append', 'delimiter', ',', 'precision', 10)



%% Output only Re and KC (if needed)
if 1
    out = [data, ReMean, ReSd, KCMean, KCSd, Re2Mean, Re2Sd, KC2Mean, KC2Sd]; 
    cd(cdata);

    % make a header of column names
    header = fgetl(fopen(infile)); fclose('all');
    header = regexprep(header, '\t', ',');
    header = [header, ',ReMean, ReSd, KCMean, KCSd,Re2Mean, Re2Sd, KC2Mean, KC2Sd'];

    % write out a csv
    dlmwrite(['VegData2_Re_', date(), '.csv'], header, 'delimiter', '', 'precision', 10)
    dlmwrite(['VegData2_Re_', date(), '.csv'], out, '-append', 'delimiter', ',', 'precision', 10)
end   
