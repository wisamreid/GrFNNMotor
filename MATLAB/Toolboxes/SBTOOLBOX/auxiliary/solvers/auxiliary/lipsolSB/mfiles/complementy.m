function [Rxz,Rsw,dgap] = complementy(x,z,s,w)
% COMPLEMENTY - Evaluate the complementarity vectors and gap.
%	Usage:  [Rxz,Rsw,dgap] = complementy(x,z,s,w)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Ubounds_exist

Rxz = x.*z; Rsw = [];
if (Ubounds_exist) Rsw = s.*w; end
dgap = full(sum([Rxz; Rsw]));
