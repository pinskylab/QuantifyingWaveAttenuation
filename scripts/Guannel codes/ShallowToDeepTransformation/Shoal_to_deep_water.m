g=9.81;rho=1024;

h=2;T=5;sig=2*pi/T;H=.29; % h and T are for Moller outer station; use h=0.5 and T=5 for inner station
k=func_disp(T,h);

%Water depth charact.
if k*h>pi
    disp(['Deep Water'])
elseif k*h<pi/10
    disp(['Shallow water'])
else disp(['Interm Water'])
end

%Energy flux
c=sig/k;
n=1/2*(1+(2*k*h(1))/sinh(2*k*h));
Cg=c*n;
E=1/8*rho*g*H^2;
Ef=E*Cg;   %En Flux

%Deep water char.
Co=g*T/(2*pi);
Cgo=1/2*Co;

%% Shoal back wave to deep water
Ho=H.*sqrt(Cg./Cgo);