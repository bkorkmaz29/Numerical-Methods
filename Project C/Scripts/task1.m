clear;clc;
format long g;

% Main for Task 1

x = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
y = [9.6459 2.5625 0.6829 0.3111 1.5471 1.1324 -0.0736 -3.7244 -12.2339 -23.4222 -43.5782];

plot(x,y,'k*')

fprintf("first order\n\n")
A = [1 -5 ; 1 -4;1 -3;1 -2;1 -1;1 0;1 1;1 2;1 3;1 4;1 5];
b = [9.6459; 2.5625; 0.6829; 0.3111; 1.5471; 1.1324; -0.0736; -3.7244; -12.2339; -23.4222; -43.5782];
At = transpose(A);
AtA = At*A;
Atb = At*b;
disp(AtA)
fprintf("condition number for 1st order Gram's")
c = cond(AtA);
disp(c);
disp(Atb)

fprintf("second order\n\n")
A = [1 -5 (-5)^2 ; 1 -4 (-4)^2;1 -3 (-3)^2;1 -2 (-2)^2;1 -1 (-1)^2;1 0 0;1 1 (1)^2;1 2 (2)^2;1 3 (3)^2;1 4 (4)^2;1 5 (5)^2];
b = [9.6459; 2.5625; 0.6829; 0.3111; 1.5471; 1.1324; -0.0736; -3.7244; -12.2339; -23.4222; -43.5782];
At = transpose(A);
AtA = At*A;
Atb = At*b;
disp(AtA)
fprintf("condition number for 2nd order Gram's")
c = cond(AtA);
disp(c);
disp(Atb)

fprintf("third order\n\n")
A = [1 -5 (-5)^2 (-5)^3 ; 1 -4 (-4)^2 (-4)^3;1 -3 (-3)^2 (-3)^3;1 -2 (-2)^2 (-2)^3;1 -1 (-1)^2 (-1)^3;1 0 0 0;1 1 (1)^2 (1)^3;1 2 (2)^2 (2)^3;1 3 (3)^2 (3)^3;1 4 (4)^2 (4)^3;1 5 (5)^2 (5)^3];
b = [9.6459; 2.5625; 0.6829; 0.3111; 1.5471; 1.1324; -0.0736; -3.7244; -12.2339; -23.4222; -43.5782];
At = transpose(A);
AtA = At*A;
Atb = At*b;
disp(AtA)
fprintf("condition number for 3rd order Gram's")
c = cond(AtA);
disp(c);
disp(Atb)

fprintf("fourth order\n\n")
A = [1 -5 (-5)^2 (-5)^3 (-5)^4 ; 1 -4 (-4)^2 (-4)^3 (-4)^3 ;1 -3 (-3)^2 (-3)^3 (-3)^4;1 -2 (-2)^2 (-2)^3 (-2)^4;1 -1 (-1)^2 (-1)^3 (-1)^4;1 0 0 0 0 ;1 1 (1)^2 (1)^3 (1)^4;1 2 (2)^2 (2)^3 (2)^4;1 3 (3)^2 (3)^3 (3)^4;1 4 (4)^2 (4)^3 (4)^4;1 5 (5)^2 (5)^3 (5)^4];
b = [9.6459; 2.5625; 0.6829; 0.3111; 1.5471; 1.1324; -0.0736; -3.7244; -12.2339; -23.4222; -43.5782];
At = transpose(A);
AtA = At*A;
Atb = At*b;
disp(AtA)
fprintf("condition number for 4th order Gram's")
c = cond(AtA);
disp(c);
disp(Atb)

fprintf("first order solution with QR method\n\n")
A = [11 0;0 110];
b = [-67.1504;-418.5014];
[Q,R] = qr(A);
d1 = R\Q.'*b;
disp(d1)

fprintf("second order solution with QR method\n\n")

A = [11 0 110 ; 0 110 0; 110 0 1958];
b = [-67.1504;-418.5014;-1298.2014];
[Q,R] = qr(A);
d2 = R\Q.'*b;
disp(d2)


fprintf("third order solution with QR method\n\n")

A = [11 0 110 0; 0 110 0 1958;110 0 1958 0;0 1958 0 41030];
b = [-67.1504;-418.5014;-1298.2014;-8698.6916];
[Q,R] = qr(A);
d3 = R\Q.'*b;
disp(d3)

fprintf("fourth order solution with QR method\n\n")
A = [11 0 110 0 1638;0 110 0 1958 1280;110 0 1958 0 35910; 0 1958 0 41030 20480;1638 1280 35910 20480 864518];
b = [-67.1504;-418.5014;-1298.2014;-8698.6916;-28356.541];
[Q,R] = qr(A);
d4 = R\Q.'*b;
disp(d4)

hold on;
x = -5:5;
y = d1(2)*x+d1(1);
plot(x,y,'b')

hold on;
x = -5:5;
y = d2(3)*x.^2+d2(2)*x+d2(1);
plot(x,y,'g')

hold on;
x = -5:5;
y = d3(4)*x.^3+d3(3)*x.^2+d3(2)*x+d3(1);
plot(x,y,'r')

hold on;
x = -5:5;
y = d4(5)*x.^4+d4(4)*x.^3+d4(3)*x.^2+d4(2)*x+d4(1);
plot(x,y,'y')
legend('point','first order','second order','third order','fourth order')

e1 = EucledianNorm(d1);
e2 = EucledianNorm(d2);
e3 = EucledianNorm(d3);
e4 = EucledianNorm(d4);

fprintf("Error of 1st order:")
disp(e1)
fprintf("Error of 2nd order:")
disp(e2)
fprintf("Error of 3rd order:")
disp(e3)
fprintf("Error of 4th order:")
disp(e4)
