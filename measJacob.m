function H = measJacob(x)
H = zeros(1,numel(x));
H(1,1) = 1+x(4);
H(1,2) = 0;
H(1,3) = -x(4);
H(1,4) = x(1)-x(3);
end
