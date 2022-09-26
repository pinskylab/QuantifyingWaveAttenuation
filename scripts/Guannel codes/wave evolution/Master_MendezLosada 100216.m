clear all;clc;
cd C:\NatCap\Codes

%% Constants
global g rho gam
g=9.81;rho=1024;

%% Bathymetry
dx=1;x=0:dx:2000-dx;lx=length(x); % Define distance vector x (every meter)

m =1/40; %Bed slope
ho=50; %Depth at offshore edge
h=ho-m.*x'; %Bathy profile
Le=length(h);

%Vegetation field
global N bv hv Cd
N=h*0;bv=h*0;hv=h*0;Cd=h*0; % initialize density, x-sectional area, height, and drag over the transect
xo=800;xs=1200; %Offshore and nearshore locations (from offshore) of vegetation
% xo=0;xs=0;%Offshore and nearshore locations (from offshore) of vegetation
[temp,x1]=min(abs(x-xo));
[temp,x2]=min(abs(x-xs));
N(x1:x2)=1200; %Shoot density  ./m^2
bv(x1:x2)=0.025; %Cross-sectional Area/m
hv(x1:x2)=0.2; %Shoot height
Cd(x1:x2)=0.1; %Drag coef

%% Wave info at offshore boundary
T=13;fp=1/T;sig=2*pi/T; % Wave period, frequency etc.
Ango=20; %Wave angle of approach
Ho=7.6/sqrt(2); %Deep water rms wave height

gam=0.42; % Wave breaking index
[H,Eta,Ef]=ThGu_MendezLosada(x,h,T,Ho,Ango);
% [H,Eta,Ef]=ThGu_Chen(x,h,T,Ho,Ango);

figure(1);col='g';
subplot 211;plot(x,H,'color',col);hold on;grid on;ylabel('Wave Height [m]','fontsize',12,'fontweight','bold')
subplot 212;plot(x,Eta,'color',col);hold on;grid on;ylabel('MWL [m]','fontsize',12,'fontweight','bold');
xlabel('X Dist [m]','fontsize',12,'fontweight','bold')

%% Estimate Runup
%Runup due to IG waves
Lo=g*T^2/(2*pi);betaf=0.0303;
Rig=sqrt(Ho*sqrt(2)*Lo*(0.563*betaf^2+0.004))/2;
%Runup due MWL
Rmwl=Eta(end);
%Total runup
Runup=1.1*(Rmwl+Rig)

