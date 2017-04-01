function [x_p,sysvar] = sysdyn(x,t,dt,iter,syspar,sysvar)
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
Kret = syspar.hydro.Kret;
Chs = syspar.hydro.Chs;
Cmr = syspar.external.Cmr;

% Evaluate Cummins integral
ktime = sysvar.k;
if iter == 1    
    thist = sysvar.t(1:ktime);
    nuhist = sysvar.nu(:,1:ktime);
    [mu] = cumminsint(Kret,thist,nuhist);
    sysvar.mu(:,ktime) = mu;
else
    mu = sysvar.mu(:,ktime);
end
% mu = [0;0];

% External loads
Fext = [0;0]; % Not implemented - keep zero for now

% System dynamics
eta_p = nu;
nu_p = (Mrb+Ainf)\(-Dl*nu-Dq*abs(nu).*nu-mu-(Chs+Cmr)*eta+Fext);
% nu_p(1) = (Mrb(1,1)+Ainf(1,1))\(-Dl(1,1)*nu(1)-mu(1)-(Chs(1,1)+Cmr(1,1))*eta(1)+Fext(1));
% nu_p(2) = 0;


%% Define y_p
x_p = [eta_p;nu_p];

