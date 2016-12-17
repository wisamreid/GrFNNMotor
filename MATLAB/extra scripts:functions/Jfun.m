function J = Jfun(x)
    J(1,1) = 1 + 3.0*cos(x(1));
    J(1,2) = -3*x(2)^2;
    J(2,1) = 4.0*x(1)^3;
    J(2,2) = 1-5.0*x(2)^4;
    return;
end