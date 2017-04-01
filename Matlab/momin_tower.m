Re_b = 3.25; % External radius, basis [m]
Re_t = 1.935; % External radius, top [m]
Ri_b = Re_b - 0.027; % Internal radius, basis [m]
Ri_t = Re_t - 0.019; % Internal radius, top [m]
L = 77.6; % Tower length [m]
rho = 8500; % "Effective" steel density [kg/m^3]

a = Re_b;
b = (Re_b-Re_t)/L;
c = Ri_b;
d = (Ri_b-Ri_t)/L;

alpha = a^2 - c^2;
beta = 2*(c*d-a*b);
gamma = b^2 - d^2;

m = rho*pi*(alpha*L+beta*L^2/2+gamma*L^3/3);

zg = (1/m)*rho*pi*(alpha*L^2/2+beta*L^3/3+gamma*L^4/4);

I = alpha*rho*pi*L^3/3 + beta*rho*pi*L^4/4 + gamma*rho*pi*L^5/5;


