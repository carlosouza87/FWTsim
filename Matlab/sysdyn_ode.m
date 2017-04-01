function [x_p] = sysdyn_ode(t,x)
global syspar iter
% x - current system state
% t - current instant of time
% dt - time step
% iter - iteration during time-step integration (e.g. iter = 1, 2, 3, 4 for 4th-order Runge-Kutta)
% syspar - structure with system parameters

%% Read state vector
eta = x(1:2);
nu = x(3:4);

%% Equations of motions
Mrb = syspar.strprop.Mrb;
Ainf = syspar.hydro.Ainf;
Dl = syspar.hydro.Dl;
Dq = syspar.hydro.Dq;
K11 = syspar.hydro.K11;
K15 = syspar.hydro.K15;
K51 = syspar.hydro.K51;
K55 = syspar.hydro.K55;
Chs = syspar.hydro.Chs;
Cmr = syspar.external.Cmr;

% Evaluate memory effects
Krad = [0;0]; % Not implemented - keep zero for now

% External loads
Fext = [0;0]; % Not implemented - keep zero for now

% System dynamics
eta_p = nu;
nu_p = (Mrb+Ainf)\(-Dl*nu-Dq*abs(nu).*nu-Krad-(Chs+Cmr)*eta+Fext);


%% Define y_p
x_p = [eta_p;nu_p];

