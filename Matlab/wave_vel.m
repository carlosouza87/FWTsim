clear all;close all;clc

T = 8;
zeta_a = 1.0;
h = 350;
g = 9.8;

omg = 2*pi/T;
eps = 1.0;
eps_tol = 0.01;
k = omg^2/g*0.01;
while abs(eps) >= eps_tol
    k = omg^2/(g*tanh(k*h));
    eps = omg^2/g - k*tanh(k*h);
end

% t = 0;
% x = -pi/2/k;
z = 0:-1:-120;

for  k1 = 1:length(z)
    % u = omg*zeta_a*cosh(k*(z+h))/sinh(k*h)*sin(omg*t-k*x);
    u(k1) = omg*zeta_a*cosh(k*(z(k1)+h))/sinh(k*h);
end

plot(u,z)
grid on
