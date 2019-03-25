function [dDV_dz]=DEdef(z, DV);
%%Dependent variables
nN=DV(1);
Tf=DV(2);
Tg=DV(3);
P=DV(4);
%% Set global variables 
global Min
global nt
global ID
global i
global OD
global Pt
global Db
global SD
global DVsol
global Lsol
global conversion
global Ai
global circ
global A
global nN0
global nA0
global nH0
global Fin
global R
global ks
global S1
global S2
global W
global kg
global Dp
%% Pressure/conversion calculations
PN=P*nN/(2.598*nN0+2*nN); %Pa
PH=3.9*PN; %Pa
PA=P*(2.23*nN0-2*nN)/(2.598*nN0+2*nN) %Pa
xN=(nN0-nN)/nN0; %conversion of N2


%u=600; %W/(m^2 K)

%%heat capacity sequence
Cpf=3000; %J/kgK
Cpg=3000; %J/kgK
deltaH=2*4.184*(-(0.54526+846.09/Tg+(459.734*10^6)/(Tg^3))*P/(101325)-5.34685*Tg-0.2525*10^-3 * Tg^2 + 1.69197*10^-6*Tg^3 - 9157.09) %Calculation of heat of reaction (per mole of ammonia therefore multilpied by 2) doi 10.1016/j.cherd.2017.10.021 
%deltaH=-92000; %J/mol N2


%% Rate and kinetics sequence
f1=1.78954*10^4 *1000/(101325^1.5*60^2) ;
f2=2.5714*10^16 * (1000/(1/101325)^0.5)/(60^2);
K1=f1*exp(-87027/(R*Tg)); %J/mol
K2=f2*exp(-198322/(R*Tg)); %J/mol
f=1*(-4.6757259 + 0.02354872*Tg + 4.687353*xN-3.463308*10^-5 *Tg^2-11.28031*xN^2+1.540881*10^-8*Tg^3+ 10.46627*xN^3) %Effectiveness factor from doi 10.1016/j.cherd.2017.10.021 (should be less than one. Use this as a check for validity)
%f=0.245 %minimum f
Rate=f*4*(K1*((PN*(PH^1.5))/PA)-K2*PA/(PH^1.5));
%% Altermative rate equation
%a=0.5;
%k=8.849*10^14*exp(-170561/(R*Tg));
%log10Ka=-2.691122*log10(Tg)-5.519265*10^-5*Tg+1.848863*10^-7*Tg^2+2001.6/Tg+26899;
%Ka=10^log10Ka
%Rate=2*k*(Ka^2*PN*(PH^3/PA^2)^a-(PA^2/PH^3)^(1-a))
%% Pressure drop sequence
eps=0.6;
a=0.1408;
b=0.03913/1000;
n=Min/nt(i);
Vf=(Min/nt(i))*R*Tg/(P)
%Vf=(b*n*P+n*R*Tg)/(3*P)+(2^(1/3)*(3*a*n^2*P-(b*n*P+n*R*Tg)^2))/(3*P*(-18*a*b*n^3*P^2-2*b^3*n^3*P^3+9*a*n^3*P*R*Tg-6*b^2*n^3*P^2*R*Tg-6*b*n^3*P*R^2*Tg^2-2*n^3*R^3*Tg^3+sqrt((-18*a*b*n^3*P^2-2*b^3*n^3*P^3+9*a*n^3*P*R*Tg-6*b^2*n^3*P^2*R*Tg-6*b*n^3*P*R^2*Tg^2-2*n^3*R^3*Tg^3)^2+4*(3*a*n^2*P-(b*n*P+n*R*Tg)^2)^3))^(1/3))-(-18*a*b*n^3*P^2-2*b^3*n^3*P^3+9*a*n^3*P*R*Tg-6*b^2*n^3*P^2*R*Tg-6*b*n^3*P*R^2*Tg^2-2*n^3*R^3*Tg^3+sqrt((-18*a*b*n^3*P^2-2*b^3*n^3*P^3+9*a*n^3*P*R*Tg-6*b^2*n^3*P^2*R*Tg-6*b*n^3*P*R^2*Tg^2-2*n^3*R^3*Tg^3)^2+4*(3*a*n^2*P-(b*n*P+n*R*Tg)^2)^3))^(1/3)/(3*2^(1/3)*P)
mu=3.332*10^-5;
rho=23; %avg
dPdL=-(1-eps)*Vf/(Dp*eps^3*Ai)*(150*(1-eps)*mu/Dp + 1.75*rho*Vf/Ai);
%%

Re=rho*   (Min*R*Tg/P)/(pi*(SD^2/4)-A)    *(i*OD*nt(i)+SD)/mu; %reynolds number feed
Reb=Dp*W/(nt(i)*mu); %reynolds number in catalyst bed (particle reynolds number)


%% Heat Transfer Coefficient Sequence
Pr=Cpg*mu/kg;
hf=0.023*Re^0.8 * Pr^0.4*(1/(OD*nt(i)+SD)/kg); %feed HTC
hp=(0.0255+0.055*W/(pi*(SD^2)/4 - nt(i)*Ai))*Tg+38; %bed to wall HTC W/m^2K
invU=1/hf + 1/hp + ID/(ks);

u=1/invU

%nA=nA0+2*nN0*xN;
%nH=3*nN;
%PA=nA/(nN+nA+nH)*P
%PN=P*nN/(nN+nA+nH)
%PH=P*nH/(nN+nA+nH)


dDV_dz=[-Rate;
    -((u*S1)/(W*Cpf))*(Tg-Tf);
    -((u*S1)/(W*Cpg))*(Tg-Tf)+(((-deltaH*S2)/(W*Cpg))*Rate);
    dPdL];
end
