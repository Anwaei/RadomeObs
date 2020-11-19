function G = dynJacob(x,u)
global V_M T
G = zeros(numel(x),numel(x));
G(1,1) = 1 + T/x(2)*V_M*cos(x(3)-x(1));
G(1,2) = T/x(2)^2*V_M*sin(x(3)-x(1));
G(1,3) = -T/x(2)*V_M*cos(x(3)-x(1));
G(1,4) = 0;
G(2,1) = -T*V_M*sin(x(3)-x(1));
G(2,2) = 1;
G(2,3) = T*V_M*sin(x(3)-x(1));
G(2,4) = 0;
G(3,1) = 0;
G(3,2) = 0;
G(3,3) = 1;
G(3,4) = 0;
G(4,1) = 0;
G(4,2) = 0;
G(4,3) = 0;
G(4,4) = 1;
end
