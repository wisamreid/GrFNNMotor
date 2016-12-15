% n = 64;  
% x0(1:n,1) = -1.9; 
% x0(2:2:n,1) = 2;
% options = optimoptions(@fsolve,'Display','iter','SpecifyObjectiveGradient',true);
% [x,F,exitflag,output,JAC] = fsolve(@bananaobj,x0,options);

function F=basicfun(x)
%     F=3.*x.^3-2*x.^2+x-7;
F = x^2;

end


