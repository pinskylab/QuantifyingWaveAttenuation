function Dif = Wave_Baldock_GamKnown2(Hmeas,fp,Cd,bv,N,hv,B,h,hd,flat,x,xflat,xflat2,pos)
global rho g Cf; up=[];

% Adjust veg height so that it isn't > depth
hvadj = min(hv, h);

%% Wave Model
sig=2*pi*fp;
dx=x(2)-x(1);lx=length(x); L=length(hd);
H=x*0;Ef=H;Db=H;Df=Db;Dv=Db;Cg=H;


H(1)=Hmeas(1);Ew=.125*rho*g*H(1)^2; 
[~,Cg(1)]=iterativek(sig,hd(1));
Ef(1)=Ew*Cg(1); %Energy flux

Lo=g/(2*pi*fp^2);Co=g/sig;Cgo=Co/2;%Deep water phase speed
hi=h(1);ki=iterativek(sig,hi);Li=2*pi/ki;%Wave # @ 1st grid pt
ni=.5.*(1+(2.*ki.*hi./sinh(2.*ki.*hi))); %to compute Cg @ 1st grid pt
Ci=Li*fp;  Cgi=Ci.*ni; %Group velocity @ 1st grid pt
Ho=H(1)*sqrt(Cgi/Cgo);%We are in intermediate water. Assume no brkg occured
So=Ho/Lo; %Transform significant wave height to rms wave height

for xx=1:lx-1 %Transform waves over sloping bottom
    Ef(xx)=0.125*rho*g*(H(xx)^2)*Cg(xx); %Ef at (xx)      
    Ef(xx+1)=Ef(xx)-dx*(Db(xx)+Df(xx)+Dv(xx)); %Ef at (xx+1) 

    k=iterativek(sig,h(xx+1));
    Cg(xx+1)=sig/k*.5.*(1+(2.*k.*h(xx+1)./sinh(2.*k.*h(xx+1)))); %Group velocity
%     Gam=0.76*k*h(xx+1)+0.29;% Ruessink Gamma
    Gam=0.5+0.4*tanh(33*So);
    H(xx+1)=sqrt(8*Ef(xx+1)/(rho*g*Cg(xx+1)));%Wave height at (xx+1)
    Hb=0.88/k*tanh(Gam*k*h(xx+1)/0.88);
    A=0.25*rho*g*fp*B;
    temp1=((Hb/H(xx+1))^3+1.5*Hb/H(xx+1))*exp(-(Hb/H(xx+1))^2);
    if isreal(H(xx+1))
        temp2=0.75*sqrt(pi)*(1-erf(Hb/H(xx+1)));nogo=0;
    else temp2=0;nogo=1;
    end
    Db(xx+1)=A*H(xx+1)^3/h(xx+1)*(temp1+temp2); %Dissipation due to brkg
    Dv(xx+1)=rho*Cd*bv(xx+1)*N(xx+1)*(k*g/(2*sig))^3/(2*sqrt(pi))*...
        (sinh(k*hvadj(xx+1))^3+3*sinh(k*hvadj(xx+1)))/...
        (3*k*cosh(k*h(xx+1))^3)*H(xx+1)^3;%Diss due to vegetation (MendezLosada2004), using corrected veg height
    % Dv2(xx+1)=rho*Cd*bv(xx+1)*N(xx+1)*(k*g/(2*sig))^3/(2*sqrt(pi))*...
    %    (sinh(k*hv(xx+1))^3+3*sinh(k*hv(xx+1)))/...
    %    (3*k*cosh(k*h(xx+1))^3)*H(xx+1)^3;%Diss due to vegetation (MendezLosada2004), using uncorrected veg height
    Df(xx+1)=rho*Cf/(16*sqrt(pi))*...
        (2*pi*fp*H(xx+1)/sinh(k*h(xx+1)))^3;%Diss due to bot friction 

    %If there is a flat portion in profile, assume no breaking occurs
    if ismember(x(xx+1),xflat2);Db(xx+1)=0;end
    if mean(flat); %if there is a flat portion
        temp=find(xflat<x(xx+1));
        if ~isempty(temp);up=temp(end);else up=1;end %loc of last flat portion
    end
    if h(xx+1)>hd(up);Db(xx+1)=0;end  % no breaking if the current water is deeper than the last flat section
end

Hpos=H(pos);if size(Hpos,1)~=size(Hmeas,1);Hpos=Hpos';end
Dif=sqrt(sum((Hpos-Hmeas).^2)/L); % compute RMS error
if nogo;Dif=999;end
  
end