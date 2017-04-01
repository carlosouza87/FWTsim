clear all;close all;clc

% Time parameters
dt = 0.125;
Time = 0:dt:2000;

% Load system properties
spar_prop;

% Initial states
eta0 = [-20;0];
% eta0 = [0;5*pi/180];
nu0 = [0;0];
x0 = [eta0;nu0];
Nsys = size(x0,1);

x = zeros(Nsys,length(Time));
x(:,1) = x0;

% Structure for online storage of simulation variables
sysvar.t = Time;
sysvar.eta = zeros(size(eta0,1),length(Time));
sysvar.eta(:,1) = eta0;
sysvar.nu = zeros(size(nu0,1),length(Time));
sysvar.nu(:,1) = nu0;
sysvar.k = 1;

% Dynamic simulation
for k1 = 2:length(Time)
    t = Time(k1);
    sysvar.k = k1-1;
    if t >= 500
        a=1;
    end
    [x(:,k1),sysvar] = RK4th(x(:,k1-1),t,dt,syspar,sysvar);
    sysvar.eta(:,k1) = x(1:2,k1);
    sysvar.nu(:,k1) = x(3:4,k1);    
end

figure(1)
plot(Time,x(1,:))
hold on
plot(Time(1:end-1),sysvar.mu(1,:)/(5e3),'r')
grid on

% figure(2)
% plot(Time,x(2,:)*180/pi)
% grid on

% [t,x] = ode45(@sysdyn_ode,[0 2000],x0);
% plot(t,x(:,1))



