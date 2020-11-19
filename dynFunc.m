function x_next = dynFunc(x,u)
global V_M T
x_next = zeros(size(x));
x_next(1) = x(1) - T/x(2)*V_M*sin(x(3)-x(1));
x_next(2) = x(2) - T*V_M*cos(x(3)-x(1));
x_next(3) = x(3) + T/V_M*u;
x_next(4) = x(4);
end
