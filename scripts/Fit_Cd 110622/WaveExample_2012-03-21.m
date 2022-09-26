%%%%% Example

%% Set up %%%%%%%%%%%
clear all;
homedir = '/Users/mpinsky/Documents/Stanford/NatCap/Coastal Protection Model/analysis/Fit_Cd 110622';
%homedir = 'P:\Users\Malin Pinsky\NatCap\Fit_Cd 110622'
cdata = homedir;
global g rho Cf; g=9.81; rho=1024; Cf = 0.003 % g=gravity; rho=water density; Cf = bottom friction

% Load data
infile='QuantforGreg_2011-07-04.txt';
cd(cdata);data=load(infile);cd(homedir);
len=length(data);tic

%% Set up and run model for one case

%kk= 364; Cd = 1.1667796; B=0.4; study='Bouma et al. 2005' % Bouma et al. 2005 field xsect
kk= 268; Cd = 2.41194267; B=0.4; study = 'Vo-Luong & Massel 2006, 2008' 

% Bathymetry
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

% Wave info
T=data(kk,4);sig=2*pi/T;fp=1/T;% Wave period, frequency etc. 
Hmeas=data(kk,8:3:29);Hmeas(isnan(Hmeas))=[]; % Wave height measured
if data(kk,3);irreg=1;else irreg=0;end %Indicates whether waves are (ir)regular. Matters to determine coeff. for transformation

% Vegetation field
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
    out = out(1:(length(out)-1));
    newx1 = pos(max(out)+1); % index into x for new start:first H measurement inside veg
    N=N(newx1:x2);bv=bv(newx1:x2);hv=hv(newx1:x2);
    Hmeas(out)=[];xd(out)=[];hd(out)=[]; % wave heights (measured)
    xd=xd-xd(1); % location of wave height measurements (in meters)
    h=h(newx1:x2);x=x(newx1:x2)-x(newx1); % depth and distance vectors
    lx=length(x); L=length(xd);
end;

h=smooth(h,round(lx/5)); if size(h,1)~=size(x,1);h=h';end; %Smooth h vector
pos=xd*0;for xx=1:L;[~,pos(xx)]=closest(x,xd(xx));end; %ID loc of H meast


% ID Location in depth vector where bed is flat and sloping
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


%% Transform wave

H = Wave_Baldock_GamKnown_H2(Hmeas(1),fp,Cd,bv,N,hv,B,h,hd,flat,x,xflat,xflat2);
H
%% Plot
%figure
plot(x, H, 'k')
ylabel('Wave height (m)')
xlabel('Distance (m)')
axis([-1, max(x)+1, min(-h), max(Hmeas)*1.3])
title(study);
set(gca, 'Color', 'w', 'TickLength', [-0.01, 0.025]) % set plotting bg to white
set(gcf, 'Color', 'w') % set figure background to white
hold on
plot(xd,Hmeas, '.k')
plot(x, -h, 'k--', 'LineWidth', 2) % depth
plot(x(hv>0), [-h(hv>0)+hv(hv>0)], 'g')
%H(end)
hold off
