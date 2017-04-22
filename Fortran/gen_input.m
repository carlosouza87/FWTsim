%% Time parameters
dt = 0.1;    % Time step [s]
ti = 0;      % Initial time [s]
tf = 2000;   % Final time [s]

%% Initial states - automatically defines number of DOFs
eta0 = [-20;0];
nu0 = [0;0];

%% Floater properties
spar_prop

%% Rotor properties
rotor_prop 

%% Genrate input files
write_sysprop
write_simpar
write_rotorprop