%% Modell för reaktor
% A=CO, B=H2, P=Metanol; Cu/ZnO/Al2O3 catalys 
% Reaktion: A+2B --> P

% Reaktor 1
clear
clc
P = 50; %atm    
T1 = 250+273; Ea = 123.5; % kJ/mol

F_tot1a= 95.53952943; F_A1a = 29.0983828; F_B1a = 58.19660071; 
F_V1 = 9.699391172 ; F_P1a = 0; % mol/s vid 250C

CpA = 30; CpB = 29.2; CpV = 35.6; CpP = 61.2; % [J/molK]

deltaH= -91.588;  % Reaktionsentalpi [kJ/mol]

k = 5.24*10^-5; % mol/kg_cat*min*atm


Ke = (3.567*10^-12)*exp(90.13/8.314*T1); %atm^-2

P_A1 = 4.7369; P_B1 = 9.4735; P_P1 = 26.8422; % Partialtryck komponenter [atm]

n = -1.3;   r = k*(P_A1*(P_B1)^2-(P_P1/Ke))*(P_A1*(P_B1)^0.5)^n;

p_W = 1120; %kg/m^3


V1=33 ;  %m^3, utifrån att vi har F_tot0= ca 95.5 mol/s

W1 = V1*p_W ;  % Katalys massa [kg]

%L/D == 2; 

W_start1 = 0 ;
W_final1 = V1*p_W;

[W1,U1]=ode45(@ode_eq,[W_start1 W_final1],[F_A1a F_B1a F_P1a T1],[],k,n,V1,Ke,P_A1,P_B1,P_P1,CpA,CpB,CpV,deltaH,F_A1a,F_B1a,F_P1a,F_V1a,T1)


F_A1b = U1(:,1);  
F_B1b = U1(:,2);
F_P1b = U1(:,3);

Tb = U1(:,4);  

XA1 = (F_A1-F_A)./F_A1;

F_A1c = U1(end,1); F_B1c = U1(end,2); F_P1c = U1(end,3);

T1c = U1(end,4);

F_tot1c = F_A1c + F_B1c + F_P1c + F_V1;

procentP1c = F_P1c/(F_tot1c);


%plot total omsättning mot katalys massa
figure(1)
plot(W1,XA1)
xlabel('W [kg]')
ylabel('Conversion of A')

%plot flöden mot katalys massa
figure(2)
plot(W1,F_A1b,W1,F_B1b,W1,F_P1b)
xlabel('W [kg]')
ylabel('Flöde [mol/s]')

%plot temp mot katalys massa
figure(3)
plot(W1,U1(:,4))
xlabel('W [kg]')
ylabel('T [K]')


% Reaktor 2
clear
clc

Ke = (3.567*10^-12)*exp(90.13/8.314*T1); %atm^-2

P_A1 = 4.7369; P_B1 = 9.4735; P_P1 = 26.8422; % Partialtryck komponenter [atm]

n = -1.3;   r = k*(P_A1*(P_B1)^2-(P_P1/Ke))*(P_A1*(P_B1)^0.5)^n;

p_W = 1120; %kg/m^3


V1=33 ;  %m^3, utifrån att vi har F_tot0= ca 95.5 mol/s

W1 = V1*p_W ;  % Katalys massa [kg]

%L/D == 2; 

W_start1 = 0 ;
W_final1 = V1*p_W;

[W1,U1]=ode45(@ode_eq,[W_start1 W_final1],[F_A1 F_B1 F_P1 T1],[],k,n,V1,Ke,P_A1,P_B1,P_P1,CpA,CpB,CpV,deltaH,F_A1,F_B1,F_P1,F_V1,T1)


F_A = U1(:,1);  
F_B = U1(:,2);
F_P = U1(:,3);

T = U1(:,4);  

XA1 = (F_A1-F_A)./F_A1;

F_A2 = U1(end,1); F_B2 = U1(end,2); F_P2 = U1(end,3);

T2 = U1(end,4);

F_tot2 = F_A2 + F_B2 + F_P2 + F_V1;

procentP2 = F_P2/(F_tot2);
