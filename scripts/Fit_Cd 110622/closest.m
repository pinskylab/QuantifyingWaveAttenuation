function [val,loc]=closest(X,Numb)

%[val,loc]=closest(X,Numb)
%This function estimates the value and position of the number in vector X
%that is the closest to the user defined value Numb.
%
%INPUTS: X: Known vector of numbers
%        Numb: User-specified number of interest
%
%OUTPUTS: val: value of number in vector X that is closest to Numb
%         loc: location of val in vector X

val=Numb*0;loc=val;
for kk=1:length(val)
    [~,loc(kk)]=min(abs(X-Numb(kk)));
end
val=X(loc);