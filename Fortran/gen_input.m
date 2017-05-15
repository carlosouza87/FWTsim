%% Time parameters
dt = 0.1;    % Time step [s]
ti = 0;      % Initial time [s]
tf = 2000;   % Final time [s]
t_clutch = 100; % Time limit for clutching spar motions [s]
t_ramp = 200; % Time limit for unitary ramp multipliying rotor thrust action on spar [s]

%% Initial states - automatically defines number of rigid-body DOFs
eta0 = [0;0]; % Initial positions
nu0 = [0;0]; % Initial velocities
Omg_rt0 = 0.01; % Initial rotor speed [rad/s] - must be greater than 0!

%% Wind velocity
Uwind = 12.0;

%% Floater properties
spar_prop

%% Rotor properties
rotor_prop 

%% Generator properties
generator_prop

%% Genrate input files
write_inpfiles
% write_sysprop
% write_simpar
% write_rotorprop
% write_generatorprop