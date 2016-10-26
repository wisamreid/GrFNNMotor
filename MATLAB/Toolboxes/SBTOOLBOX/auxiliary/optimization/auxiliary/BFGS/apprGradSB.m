function [grad] = apprGradSB(x,FUN,TOL)
% Approximation of the gradient of the function FUN at the point x with
% error o(TOL^2))
%
% USAGE:
% ~~~~~~ 
% [grad] = apprGradSBSBAO(x,FUN,TOL)
%
% Input Arguments:
% ~~~~~~~~~~~~~~~~
% x: Point at which the gradient is to be evaluated
% FUN: Function of which the gradient is to be calculated
% TOL: Error tolerance
%
% Output Arguments:
% ~~~~~~~~~~~~~~~~~
% grad: Calculated gradient as column vector

ndim = length(x);
grad = zeros(ndim,1);
for i=1:ndim
    e = zeros(ndim,1);
    e(i) = 1;
    grad(i) = (feval(FUN,x + TOL*e) - feval(FUN,x - TOL*e))/(2*TOL);
end
return