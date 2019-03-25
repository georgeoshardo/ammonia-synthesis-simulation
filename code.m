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
%% Pipe properties
nt=[500:100:1900];              %number of tubes
 % top temperature K
%Ttop=765 %max Ttop at 200 bar
ID=0.1016;                      %ID of catalyst pipe 0.105=4"
circ=pi*ID;                     %circumference of each pipe
i=1;                            %iteration initial index
%for i=1:length(nt);            %Commented out once optimisation complete
Ai=(pi*ID^2)/4;                 %individual pipe area m^2
A=Ai*nt(i);                     %m^2 flow area of catalyst 
thickness=0.0060198;            %m
OD=ID+thickness;                %Outer diameter
Pt=1.15*OD;                     %Tube pitch (triangular)
Db=2*((nt(i)*(OD^2)/4)/(((pi*OD^2)/(2*sqrt(3)*Pt^2))))^(1/2); %bundle diameter for triangular pitch
SD=Db+0.1;                      %Bundle thickness from fig12.10 pg 646 C&R Vol6  (m)
S1=circ*nt(i);                  % surface area of cooling tube per unit length m
S2=Ai*nt(i);                    %catalyst zone csa m^2
Dp=5/1000;                      %catalyst particle dia (m)


%%Flow properties
Min=(70602.3*1000/2.2)/(60^2);  %mol/s
Fin=Min/A;                      %Flow FLUX mol/m^2 s
W=(1041601/2.2)/(60^2);         %mass flow rate kg/s
nN0=Fin*0.1485434;
nA0=Fin*0.3388929;
nH0=Fin*0.3634034;

%% Phyiscal constants
R=8.314
ks=50 %W/mK (stainless steel conductivity)
kg=0.050; %W/mK thermal conductivity of gas

% Specific heat constants (kcal/kmol) from doi 10.1016/j.cherd.2017.10.021
%Hydrogen
AH=6.952;
BH=-0.04576*10^2;
CH=0.09563*10^5;
DH=-0.2079*10^9
%Nitrogen
AN=6.903;
BN=-0.03753*10^2;
CN=0.193*10^5;
DN=-0.6861*10^9;
%Methane
AM=4.750;
BM=1.2*10^2;
CM=0.303*10^5;
DM=-0.263*10^9;



%% Solving the ODE system
domain=[0 14];
Ttop=715
initialconditions=[nN0 Ttop Ttop 100*10^5];
[Lsol, DVsol]=ode45(@DEdef,domain,initialconditions);

conversion=(nN0-DVsol(:,1))/nN0;            %conversion profile
plotgraphs %Plot graphs function (see plotgraphs.m)



%%Calculating reactor properties (the interesting bit) 
vel=(Min*R*DVsol(:,3)/(1.5*10^7))/A;        %velocity
reactorlength=Lsol((find(conversion>0.2,1))) %length function
catalystvolume=reactorlength*A              %Catalyst volume 
volflow=Min*8.314*DVsol(:,3)/(175*10^5);    %volume flow profile
totalvolume=pi*(SD^2)/4 * reactorlength     %total reactor volume (feed volume + catalyst volume)
spacevel=volflow/catalystvolume*60^2;        %space velocity per hour
tau=catalystvolume./volflow;                 %residence time s


datafor=[reactorlength;
    catalystvolume;
    totalvolume]



