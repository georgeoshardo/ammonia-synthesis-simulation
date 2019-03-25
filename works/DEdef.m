function [dDV_dz]=DEdef(z, DV);

nN=DV(1);
Tf=DV(2);
Tg=DV(3);
P=DV(4);

nt=900 %number of tubes
ID=0.15   %ID of catalyst pipe m
Ai=(pi*ID^2)/4 %individual pipe area m^2
circ=pi*ID

Min=70602.3/2.2; %kmol/hr
A=Ai*nt; %m^2 flow area of catalyst 
Fin=Min/A;

nN0=Fin*0.1485434; %kmol/m2hr
nA0=Fin*0.3388929; %kmol/m2hr
nH0=Fin*0.3634034; %kmol/m2hr

%P=130; %Pa
R= 1.987; %kcal/K mol




PN=P*nN/(2.598*nN0+2*nN); %bar
PH=3*PN; %bar
PA=P*(2.23*nN0-2*nN)/(2.598*nN0+2*nN); %bar

f=1;
f1=1.78954*10^4;
f2=2.5714*10^16;

K1=f1*exp(-20800.2/(R*Tg)); %kcal/kmol
K2=f2*exp(-47400/(R*Tg)); %kcal/kmol


u=900; %kcal/(hr m^2 K)
S1=circ*nt; %m^2
W=1041601/2.2; %kg/hr
Cpf=0.79; %kcal/kgK
Cpg=0.707; %kcal/kgK
deltaH=-26600; %kcal/kmol N2
S2=Ai*nt; %m^2

Rate=f*(K1*((PN*(PH^1.5))/PA)-K2*PA/(PH^1.5));

eps=0.3;
Vf=Min*(8.3144598*10^-8)*Tg/(P*60^2);
Dp=5/1000;
Ac=A;
mu=3.81*10^-5;
rho=0.5;

dPdL=-(1-eps)*Vf/(Dp*eps^3*A)*(150*(1-eps)*mu/Dp + 1.75*rho*Vf/A)

dDV_dz=[-Rate;
    -((u*S1)/(W*Cpf))*(Tg-Tf);
    -((u*S1)/(W*Cpg))*(Tg-Tf)+(((-deltaH*S2)/(W*Cpg))*Rate);
    dPdL];
end