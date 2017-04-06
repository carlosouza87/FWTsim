%% Rigid body inertia matrix
m = 7466330 + 249718 + 24000 + 54000; % hull + tower + nacelle + hub [kg]
zg_spar = -89.92; % z coordinate for spar's center of gravity [m]
zg_tower = 43.4; % z coordinate for tower's center of gravity [m]
zg_nac = 87.6; % z coordinate for nacelle's center of gravity [m]
zg_hub = 87.6; % z coordinate for hub's center of gravity [m]
Iyy = 6.46e10 + 4.22e8 + 4.14e8 + 1.84e8; % hull + tower + nacelle + hub [kg.m^2]
zg = (zg_spar*7466330+zg_tower*249718+zg_nac*24000+zg_hub*54000)/m; % [m]

% Mrb = [m m*zg;m*zg Iyy]; % Rigid-body inertia matrix
Mrb = [m 0;0 Iyy]; % Decoupled system
syspar.strprop.Mrb = Mrb;

%% Infinite-frequency added mass matrix;
% Ainf = [7.76e6 -4.83e8;-4.83e8 3.79e8];
% Ainf = [7.76e6 0;0 3.79e8]; % Decoupled system
Ainf = [7.99e6 0;0 3.81e10]; % Decoupled system, and with A0 instead
syspar.hydro.Ainf = Ainf;

%% Linear damping matrix
Dl = [100000 0;0 0];
syspar.hydro.Dl = Dl;

%% Quadratic damping matrix (Morison viscous term)
D1 = 6.50; % Diameter of spar's upper submerged part [m]
D2 = 7.95; % (Mean) diameter of spar's intermediate (tappered) submerged part [m]
D3 = 9.40; % Diameter of spar's lower submerged part [m]

L1 = 4.0; % Length of spar's upper submerged part [m]
L2 = 8.0; % Length of spar's intermediate (tappered) submerged part [m]
L3 = 108.0; % Length of spar's lower submerged part [m]

h1 = 2.0; % Depth of center of spar's upper submerged part [m]
h2 = 8.0; % Depth of center of spar's intermediate submerged part [m]
h3 = 108.0; % Depth of center of spar's lower submerged part [m]

Cd = 0.6; % Morison viscous coefficient

% Quadratic damping coefficients WITHOUT REGARD OF WAVE PARTICLE VELOCITY!
d11q = 0.5*Cd*(D1*L1+D2*L2+D3*L3);
d55q = 0.5*Cd*(D1*L1*h1+D2*L2*h2+D3*L3*h3);

Dq = [d11q 0;0 d55q]; % Quadratic damping matrix
syspar.hydro.Dq = Dq;

%% Retardation functions
% load mem_11.txt
% load mem_15.txt
% load mem_55.txt
% 
% syspar.hydro.Kret(1,1).T = mem_11(:,1);
% syspar.hydro.Kret(1,1).K = mem_11(:,2);
% syspar.hydro.Kret(1,2).T = mem_15(:,1);
% syspar.hydro.Kret(1,2).K = mem_15(:,2);
% syspar.hydro.Kret(2,1).T = mem_15(:,1);
% syspar.hydro.Kret(2,1).K = mem_15(:,2);
% syspar.hydro.Kret(2,2).T = mem_55(:,1);
% syspar.hydro.Kret(2,2).K = mem_55(:,2);

% dTret = 0.25;

%% Hydrostatic stiffness matrix
Chs = [0 0;0 1.59e9];

syspar.hydro.Chs = Chs;

%% Mooring stiffness matrix
% Cmr = [41180 -2816000;-2816000 311100000];
Cmr = [41180 0;0 311100000];  % Decoupled system

syspar.external.Cmr = Cmr;


