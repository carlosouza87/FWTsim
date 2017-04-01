function [f,S] = pspec(t,x)

L = length(t);
dt = t(2) - t(1);
df = 1/dt;
X = fft(x);

S2 = abs(X/L);
S = S2(1:round(L/2)+1);
S(2:end-1) = 2*S(2:end-1);

f = df*(0:round(L/2))/L;