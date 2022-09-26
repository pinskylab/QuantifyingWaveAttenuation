%This function computes wave number based on water depth and wave period
%using the linear dispersion relation and an iterative Newton-Raphson
%technique

function [k2]=func_disp(T,h);

%Define Constants
sigma=2*pi/T;
g=9.81;
k=sigma^2/g*(tanh(sigma^2*h/g)).^(-0.5);%rough approx of k
X(1)=k*h;

cst=sigma^2*h/g;

%Initialize loop
i=1;
f=cst-X(1)*tanh(X(1));
r=-tanh(X(1))-X(1)*(1-tanh(X(1))^2);

%Loop to determine f

while abs(f)>=1e-3
    i=i+1;
    X(i)=X(1)-f/r;
    f=cst-X(i)*tanh(X(i));
    r=-tanh(X(i))-X(i)*(1-tanh(X(i))^2);
end;

k2=X(i)/h;

    
