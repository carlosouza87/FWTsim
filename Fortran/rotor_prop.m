%% Blade properties
% Names of files with foils data (order should be in accordance with blade geometry!)
foilfile(1,1) = {'cyl1.dat'};
foilfile(2,1) = {'cyl2.dat'};
foilfile(3,1) = {'DU40_A17.dat'};
foilfile(4,1) = {'DU35_A17.dat'};
foilfile(5,1) = {'DU30_A17.dat'};
foilfile(6,1) = {'DU25_A17.dat'};
foilfile(7,1) = {'DU21_A17.dat'};
foilfile(8,1) = {'NACA64_A17.dat'};
Bl = 3; % Number of blades []
Irt = 115.926E3; % Blades + hub inertia about rotor axis [kg.m^2]
Rtip = 63.0; % Rotor radius [m]
Rhub = 1.5; % Hub radius [m]
Zhub = zg_hub; % Hub height [m]
% Blade dimensions:
% Local radius [m] / Twist angle [deg] / Element length [m] / Chord [m] / Foil id []
Blade_dim = [2.8667 13.308 2.7333 3.542 1; 
    5.6000 13.308 2.7333 3.854 1; 
    8.3333 13.308 2.7333 4.167 2; 
    11.7500 13.308 4.1000 4.557 3;
    15.8500 11.480 4.1000 4.652 4;
    19.9500 10.162 4.1000 4.458 4;
    24.0500 9.011 4.1000 4.249 5; 
    28.1500 7.795 4.1000 4.007 6; 
    32.2500 6.544 4.1000 3.748 6; 
    36.3500 5.361 4.1000 3.502 7; 
    40.4500 4.188 4.1000 3.256 7; 
    44.5500 3.125 4.1000 3.010 8; 
    48.6500 2.319 4.1000 2.764 8; 
    52.7500 1.526 4.1000 2.518 8; 
    56.1667 0.863 2.7333 2.313 8; 
    58.9000 0.370 2.7333 2.086 8; 
    61.6333 0.106 2.7333 1.419 8]; 