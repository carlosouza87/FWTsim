close all
load results.txt -ascii

Time = results(:,1);
lt = length(Time);
% k_13 = round(lt/3);
kclutch = 100/0.1+1;
x1 = results(:,2);
[fx1,Sx1] = pspec(Time(kclutch:end),x1(kclutch:end)-mean(x1(kclutch:end)));
x5 = results(:,3);
[fx5,Sx5] = pspec(Time(kclutch:end),x5(kclutch:end)-mean(x5(kclutch:end)));
x1d = results(:,4);
[fx1d,Sx1d] = pspec(Time(kclutch:end),x1d(kclutch:end)-mean(x1d(kclutch:end)));
x5d = results(:,5);
[fx5d,Sx5d] = pspec(Time(kclutch:end),x5d(kclutch:end)-mean(x5d(kclutch:end)));
Omg_rt = results(:,6);
[fOmg_rt,SOmg_rt] = pspec(Time(kclutch:end),Omg_rt(kclutch:end)-mean(Omg_rt(kclutch:end)));
beta = results(:,7);
[fbeta,Sbeta] = pspec(Time(kclutch:end),beta(kclutch:end)-mean(beta(kclutch:end)));
Th = results(:,8);
[fTh,STh] = pspec(Time(kclutch:end),Th(kclutch:end)-mean(Th(kclutch:end)));
Qaero = results(:,9);
[fQaero,SQaero] = pspec(Time(kclutch:end),Qaero(kclutch:end)-mean(Qaero(kclutch:end)));
Qgen = results(:,10);
[fQgen,SQgen] = pspec(Time(kclutch:end),Qgen(kclutch:end)-mean(Qgen(kclutch:end)));

figure(1)
% subplot(1,2,1)
subplot(2,1,1)
plot(Time,x1)
xlabel('Time (s)')
ylabel('x1 (m)')
title('FWT motion - surge')
grid on
% subplot(1,2,2)
subplot(2,1,2)
plot(1./fx1,Sx1)
xlabel('Period (s)')
ylabel('S_{x1} (m^2.s)')
title('PSD - surge')
xlim([0 200])
grid on

figure(2)
subplot(2,1,1)
plot(Time,x5*180/pi)
xlabel('Time (s)')
ylabel('x5 (deg)')
title('FWT motion - pitch')
grid on
subplot(2,1,2)
plot(1./fx5,Sx5)
xlabel('Period (s)')
ylabel('S_{x5} (rad^2.s^-1)')
title('PSD - pitch')
xlim([0 200])
grid on

figure(3)
subplot(2,1,1)
plot(Time,x1d)
xlabel('Time (s)')
ylabel('x1p (m/s)')
title('FWT velocity - surge')
grid on
subplot(2,1,2)
plot(1./fx1d,Sx1d)
xlabel('Period (s)')
ylabel('S_{x1} (m^2.s-1)')
title('PSD - surge velocity')
xlim([0 200])
grid on

figure(4)
subplot(2,1,1)
plot(Time,x5d*180/pi)
xlabel('Time (s)')
ylabel('x5p (deg/s)')
title('FWT velocity - pitch')
grid on
subplot(2,1,2)
plot(1./fx5d,Sx5d)
xlabel('Period (s)')
ylabel('S_{x1} (rad^2.s^-1)')
title('PSD - pitch velocity')
xlim([0 200])
grid on

figure(5)
subplot(2,1,1)
plot(Time,Omg_rt*180/pi)
xlabel('Time (s)')
% ylabel('\omega_{rotor} (RPM)')
ylabel('\omega_{rotor} (deg/s)')
title('Rotor speed')
grid on
subplot(2,1,2)
plot(1./fOmg_rt,SOmg_rt)
xlabel('Period (s)')
ylabel('S_{Omg_{rt}} (rad^2.s^-1)')
title('PSD - rotor Speed')
xlim([0 200])
grid on

figure(6)
subplot(2,1,1)
plot(Time,beta*180/pi)
xlabel('Time (s)')
ylabel('\beta (deg)')
title('Blade pitch angle')
grid on
subplot(2,1,2)
plot(1./fbeta,Sbeta)
xlabel('Period (s)')
ylabel('S_{beta} (rad^2.s^-1)')
title('PSD - blade pitch angle')
xlim([0 200])
grid on

figure(7)
subplot(2,1,1)
plot(Time,Th)
xlabel('Time (s)')
ylabel('Thrust (N)')
title('Rotor thrust')
grid on
subplot(2,1,2)
plot(1./fTh,STh)
xlabel('Period (s)')
ylabel('S_{Th} (rad^2.s^-1)')
title('PSD - Rotor thrust')
xlim([0 200])
grid on

figure(8)
subplot(2,1,1)
plot(Time,Qaero)
xlabel('Time (s)')
ylabel('Torque (N.m)')
title('Rotor torque')
grid on
subplot(2,1,2)
plot(1./fQaero,SQaero)
xlabel('Period (s)')
ylabel('S_{Qaero} (rad^2.s^-1)')
title('PSD - rotor torque')
xlim([0 200])
grid on

% figure(9)
% subplot(2,1,1)
% plot(Time,Omg_rt*97)
% hold on
% plot(Time,121.6805*ones(1,length(Time)),'r--')
% xlabel('Time (s)')
% % ylabel('\omega_{rotor} (RPM)')
% ylabel('\omega_{gen} (rad/s)')
% title('Generator speed')
% grid on
% subplot(2,1,2)
% plot(1./fQaero,SQaero)
% xlabel('Period (s)')
% ylabel('S_{Qaero} (rad^2.s^-1)')
% title('PSD - rotor torque')
% xlim([0 200])
% grid on


