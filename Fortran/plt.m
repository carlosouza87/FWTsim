close all
load results.txt -ascii

Time = results(:,1);
x1 = results(:,2);
x5 = results(:,3);
x1d = results(:,4);
x5d = results(:,5);
Omg_rt = results(:,6);
Th = results(:,7);
Qaero = results(:,8);

figure(1)
plot(Time,x1)
xlabel('Time (s)')
ylabel('x1 (m)')
title('FWT motion - surge')
grid on

figure(2)
plot(Time,x5*180/pi)
xlabel('Time (s)')
ylabel('x5 (deg)')
title('FWT motion - pitch')
grid on

figure(3)
plot(Time,x1d)
xlabel('Time (s)')
ylabel('x1p (m)')
title('FWT motion - surge')
grid on

figure(4)
plot(Time,x5d)
xlabel('Time (s)')
ylabel('x5p (m)')
title('FWT motion - pitch')
grid on

figure(5)
plot(Time,Omg_rt*60/(2*pi))
xlabel('Time (s)')
ylabel('\omega_{rotor} (RPM)')
title('Rotor speed')
grid on

figure(6)
plot(Time,Th)
xlabel('Time (s)')
ylabel('Thrust (N)')
title('Rotor thrust')
grid on

figure(7)
plot(Time,Qaero)
xlabel('Time (s)')
ylabel('Moment (N.m)')
title('Rotor torque')
grid on
