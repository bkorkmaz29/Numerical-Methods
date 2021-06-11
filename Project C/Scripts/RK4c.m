function [x1, x2] = RK4c(f1, f2, h, iX1, iX2)
% function for Runge-Kutta method 4th order with constant step
x1(1) = iX1;
x2(1) = iX2;
for i = 1:1:(20 / h)
 [k1, k2] = Ks(f1, f2, x1(i), x2(i), h);
 x1(i + 1) = x1(i) + h * (k1(1) + 2 * k1(2) + 2 * k1(3) + k1(4)) / 6;
 x2(i + 1) = x2(i) + h * (k2(1) + 2 * k2(2) + 2 * k2(3) + k2(4)) / 6;
end
