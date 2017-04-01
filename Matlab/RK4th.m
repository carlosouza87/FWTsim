function [xk1,sysvar] = RK4th(xk,t,dt,syspar,sysvar)
% 4-th order Runge-Kutta integration

iter = 1;
[k1,sysvar] = sysdyn(xk,t,dt,iter,syspar,sysvar);
% k1 = dt*sysdyn(xk,t,dt,iter,syspar,sysvar);
iter = 2;
[k2,~] = sysdyn(xk+k1*dt/2,t+dt/2,dt,iter,syspar,sysvar);
% k2 = dt*sysdyn(xk+k1/2,t+dt/2,dt,iter,syspar,sysvar);
iter = 3;
[k3,~] = sysdyn(xk+k2*dt/2,t+dt/2,dt,iter,syspar,sysvar);
% k3 = dt*sysdyn(xk+k2/2,t+dt/2,dt,iter,syspar,sysvar);
iter = 4;
[k4,~] = sysdyn(xk+k3*dt,t+dt,dt,iter,syspar,sysvar);
% k4 = dt*sysdyn(xk+k3,t+dt,dt,iter,syspar,sysvar);

xk1 = xk + (dt/6)*(k1+2*k2+2*k3+k4);
% xk1 = xk + (1/6)*(k1+2*k2+2*k3+k4);