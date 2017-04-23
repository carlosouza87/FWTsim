clear all;close all;clc
dos('dbg')
% dos('FWTsim')
load results.txt

figure(1)
plot(results(:,1),results(:,2))
grid on

figure(2)
plot(results(:,1),results(:,6))
grid on