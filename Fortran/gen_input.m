%% Time parameters
dt = 0.1;    % Time step [s]
ti = 0;      % Initial time [s]
tf = 2500;   % Final time [s]
t_clutch = 500; % Time limit for clutching system dynamics [s]

%% Initial states - automatically defines number of DOFs
eta0 = [-20;0];
nu0 = [0;0];

%% Wind velocity
Uwind = 11.4;

%% Floater properties
spar_prop

%% Rotor properties
rotor_prop 

%% Genrate input files
write_sysprop
write_simpar
write_rotorprop