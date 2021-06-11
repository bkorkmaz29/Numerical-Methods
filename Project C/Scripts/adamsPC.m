function [x1, x2] = adamsPC(f1, f2, h, iX1, iX2)
% Function for Adams PC method
x1(1) = iX1;
x2(1) = iX2;
betaE = [1901/720, -2774/720, 2616/720, -1274/720, 251/720];
betaI = [475/1440, 1427/1440, -798/1440, 482/1440, -173/1440, 27/1440];
for i = 1:4
 [k1, k2] = Ks(f1, f2, x1(i), x2(i), h);
 x1(i + 1) = x1(i) + h * (k1(1) + 2 * k1(2) + 2 * k1(3) + k1(4)) / 6;
 x2(i + 1) = x2(i) + h * (k2(1) + 2 * k2(2) + 2 * k2(3) + k2(4)) / 6;
end
for i = 6:ceil(20 / h)
 sumX1 = 0;
 sumX2 = 0;
 for j = 1:5
 sumX1 = sumX1 + betaE(j) * f1(x1(i - j), x2(i - j));
 sumX2 = sumX2 + betaE(j) * f2(x1(i - j), x2(i - j));
 end
 tempX1 = x1(i - 1) + h * sumX1;
 tempX2 = x2(i - 1) + h * sumX2;
 sumX1 = 0;
 sumX2 = 0;
 for j = 1:5
 sumX1 = sumX1 + betaI(j + 1) * f1(x1(i - j), x2(i - j));
 sumX2 = sumX2 + betaI(j + 1) * f2(x1(i - j), x2(i - j));
 end
 x1(i) = x1(i - 1) + h * sumX1 + h * betaI(1) * f1(tempX1, tempX2);
 x2(i) = x2(i - 1) + h * sumX2 + h * betaI(1) * f2(tempX1, tempX2);
end
end