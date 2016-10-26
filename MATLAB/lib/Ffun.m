function [F] = Ffun(x)

    F(1,1) = x(1)-x(2)^3+3.0*sin(x(1));
    F(2,1) = x(2)+x(1)^4-x(2)^5;
    return;
end