function [k,Cg,n]=iterativek(sigma,hdummy,Uo)

%function [k]=iterativek(sigma,hdummy)
% The following loop performs an iterative solution for k 
% The calculation is based on a given and constant angular wave frequency
warning off MATLAB:divideByZero
if nargin < 3
	Uo = 0; end
	
g=9.81;                                              % Gravity (m/s^2)
b=((sigma.^2).*hdummy)./g;                           % (dimensionless)
kestimated=((sigma.^2)./(g.*(((tanh(b)).^(1/2)))));  % Estimated wavenumber (m^-1) to begin the loop 
                                                     % (based onEckarts approximation.
kprevious=0.0000001; 
count = 0;
while  ((abs(kestimated-kprevious) > 0.000005) & (count < 1000));  % Designates 0.000005 as acceptable error
	count = count+1;
   kh=(kestimated.*hdummy);                        % First calculates kh using the estimated k value
   cur = (sigma.^2).*(1-Uo.*kestimated./sigma).^2;
   kcalculated= cur./(tanh(kh).*g);      % Using the kh value the loop calculates k
   kprevious=kestimated;                        % The loop now sets the 'previous' k value to the 'estimated'
   kestimated=kcalculated;                      % And it sets the 'estimated' value to the 'calculated'
                                                % The loop then continues by subtracting the
                                                % two values and testing for the error between them.
end

if isnan(kestimated);kcalculated=NaN;end
k=kcalculated;

if size(k,1)~=size(hdummy,1);k=k';end

n=.5.*(1+(2.*k.*hdummy./sinh(2.*k.*hdummy)));
C=sigma/k;  Cg=C.*n; %Group velocity




