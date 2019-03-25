nt=900 %number of tubes
ID=0.15   %ID of catalyst pipe
Ai=(pi*ID^2)/4 %individual pipe area m^2
A=Ai*nt; %m^2 flow area of catalyst 


Min=70602.3/2.2; %kmol/hr
Fin=Min/A;

nN0=Fin*0.1485434;
nA0=Fin*0.3388929;
nH0=Fin*0.3634034;
Tin=770; %K

domain=[0 20];
initialconditions=[nN0 Tin Tin 150];

[Lsol, DVsol]=ode45(@DEdef,domain,initialconditions);

subplot(3,1,1)
hold on
plot(Lsol,(DVsol(:,2)))
plot(Lsol,DVsol(:,3))
legend('Feed gas','Reacting gas')
xlabel('Distance through reactor')
ylabel('Temperature (K)')
hold on
subplot(3,1,2)
plot(Lsol,(nN0-DVsol(:,1))/nN0)
xlabel('Distance through reactor')
ylabel('N_2 conversion')
subplot(3,1,3)
plot(Lsol,DVsol(:,4))
ylabel('Pressure (atm)')
xlabel('Distance through reactor')