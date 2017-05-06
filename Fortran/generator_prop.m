%% Drivetrain properties
Igen = 534.116; % Generator inertia about HSS [kg.m^2]
Ngr = 97.0; % Gear ratio []

%% Blade pitch controller properties
PC_KI = 0.0008965149; % Integral gain for pitch controller [s]
PC_KK = 0.1099965; % Pitch angle were the the derivative of the aerodynamic power w.r.t. pitch has increased by a factor of two relative to the derivative at rated pitch (zero) [rad]
PC_KP = 0.006275604; % Proportional gain for pitch controller at rated pitch (zero) [s]
PC_MaxPit = 1.570796; % Maximum pitch setting in pitch controller [rad]
PC_MaxRat = 0.1396263; % Maximum pitch  rate (in absolute value) in pitch  controller [rad/s]
PC_MinPit = 0.0; % Minimum pitch setting in pitch controller [rad]
PC_RefSpd = 122.9096; % Desired (reference) HSS speed for pitch controller [rad/s]

%% Generator torque controller properties
VS_CtInSp = 70.16224; % Transitional generator speed (HSS side) between regions 1 and 1 1/2 [rad/s]
VS_DT = 0.00125; % Communication interval for torque controller [s]
VS_MaxRat = 15000.0; % Maximum torque rate (in absolute value) in torque controller [N.m/s]
VS_MaxTq = 47402.91; % Maximum generator torque in Region 3 (HSS side) [N.m]
VS_Rgn2K = 2.332287; % Generator torque constant in Region 2 (HSS side) [N.m.s^2/rad^2]
VS_Rgn2Sp = 91.21091; % Transitional generator speed (HSS side) between regions 1 1/2 and 2 [rad/s]
VS_Rgn3MP = 0.01745329; % Minimum pitch angle at which the torque is computed as if we are in region 3 regardless of the generator speed [rad]
VS_RtGnSp = 121.6805; % Rated generator speed (HSS side) [rad/s]
VS_RtPwr = 5296610.0; % Rated generator generator power in Region 3 [W]
VS_SlPc = 10.0; % Rated generator slip percentage in Region 2 1/2 [%]