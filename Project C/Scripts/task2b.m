clc; clear all; close all;
 
f = @(x) [x(2) + x(1) * (0.5 - x(1) .^ 2 - x(2) .^ 2); -x(1) + x(2) * (0.5 - x(1) .^ 2 - x(2) .^ 2)];
fode45 = @(t,x) f(x);
 
a = 0;
b = 20;
x0 = [-0.002; -0.02];
eps_abs = 0.001;
eps_rel = 0.001;
beta = 0.9;
 
[t_ode45,y_ode45] = ode45(fode45, [a;b], x0);
 
f1 = @(t,x) x(2) + x(1) * (0.5 - x(1) .^ 2 - x(2) .^ 2);
f2 = @(t,x) -x(1) + x(2) * (0.5 - x(1) .^ 2 - x(2) .^ 2);
[t,y,err,H] = rk4d(f, x0, 0.01, 0.0000001, a, b, eps_rel, eps_abs, beta);
 
figure;
plot(y(:,1), y(:,2), y_ode45(:,1), y_ode45(:,2));
xlabel('x1'); ylabel('x2');
legend('RK4','ode45');
title('RK4 with step size automatically adjusted vs. ode45');
 
figure;
plot(t,H);
xlim([a b]);
xlabel('t'); ylabel('step size');
title('step size vs. time');
 
figure;
plot(t,err);
xlim([a b]);
xlabel('t'); ylabel('error');
title('error estimate vs. time');
 
figure;
plot(t,y);
xlim([a b]);
xlabel('t'); ylabel('x');
title('solution vs. time');
legend('x1','x2');
