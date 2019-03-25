function plotgraphs

global DVsol
global Lsol
global conversion

subplot(4,1,1)
hold on
plot(Lsol,(DVsol(:,2)))
plot(Lsol,DVsol(:,3))
legend('Feed gas','Reacting gas')
xlabel('Distance through reactor (m)')
ylabel('Temperature (K)')
ylim([500 800])
hold on
subplot(4,1,2)
plot(Lsol,conversion)
xlabel('Distance through reactor (m)')
ylabel('N_2 conversion')
subplot(4,1,3)
plot(Lsol,DVsol(:,4)/10^5)
ylabel('Pressure (Bar)')
xlabel('Distance through reactor (m)')
subplot(4,1,4)
plot(Lsol,DVsol(:,1))
ylabel('Rate (mol N_2 / m^3 s)')
xlabel('Distance through reactor (m)')
end