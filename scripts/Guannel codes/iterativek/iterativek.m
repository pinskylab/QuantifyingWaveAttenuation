function [k]=iterativek(sigma,hdummy,Uo)

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

% L=(2*pi)/kcalculated
% check=h/L








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               ACCURACY CHECK                                                     %
% If the iterative solution for k is correct then for the entire range of kh I should satisfy      %
% the dispersion relation.  Therefore sigma^2=g*k*tanh(kh).  Since the attached plot is a constant %
% line at zero I have satisfied the equation for an array of 75 h and k values (a substantial data %
% array. Therefore it is safe to say that the iterative solution is accurate.  The degree of error %
% in the program can be seen by changing the scale of the plot (see figure 2).  Although not       %
% quantified as a percent error it indicates that the degree of error is around 10^-17             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(1);
%plot(kh(1:75),((sigma.^2)-(g.*kcalculated(1:75).*tanh(kh(1:75))))),axis([0 1.5 -1 5])
%xlabel('kh');
%title('Figure 1: Dispersion Relation vs. kh (Accuracy Check)');

%figure(2);
%plot(kh(1:75),((sigma.^2)-(g.*kcalculated(1:75).*tanh(kh(1:75)))))
%xlabel('kh');
%title('Figure 2: Dispersion Relation vs. kh (Degree of Error)');
