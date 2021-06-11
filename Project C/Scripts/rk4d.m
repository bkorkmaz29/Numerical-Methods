function [t,y,e,h] = rk4d(f, x, h, h_min, a, b, eps_rel, eps_abs, beta)
% function for RK4 method with adjusted step 
Y = zeros(int32((b-a)/h_min),2); Y(1,:) = x;
T = zeros(int32((b-a)/h_min),1); T(1) = a;
E = zeros(int32((b-a)/h_min),1); E(1) = 0;
H = zeros(int32((b-a)/h_min),1); H(1) = h;
t = a;
j = 2;

while t < b
  tmp1 = rk4n(f,x,h,2);
  s1 = tmp1(2,:);
  tmp2 = rk4n(f,x,h/2,3);
  s2 = tmp2(3,:);
  delta = abs(s2 - s1) / 15.0;
  epsilon = eps_abs + abs(s2) * eps_rel;
  r = epsilon ./ (delta .^ 0.2);
  if r(1) < r(2)
    d = r(1);
    err = epsilon(1);
  else
    d = r(2);
    err = epsilon(2);
  end 
  if beta * d >= 1.0
    t = t + h;
    T(j) = t;
    H(j) = h;
    x = s1';
    Y(j,:) = x;
    E(j) = err;
    j = j + 1;
  end
  h = beta * h * d;
  if h < h_min
    error('method failed with this h_min');
    return;
  end
end 

y = Y(1:(j-1),:);
t = T(1:(j-1),:);
e = E(1:(j-1),:);
h = H(1:(j-1),:);

