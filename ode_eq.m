function dUdW = ode_func(W,U,k,V,n,Ke,P_A,P_B,P_P,CpA,CpB,CpV,deltaH,F_A,F_B,F_P,F_V,T)
    
%extract FA, FB  and FP from U
F_A = U(1);  F_B = U(2); F_P = U(3);  
T = U(4);

%calculate concentrations
CA=F_A/V;  CB=F_B/V;

deltaH = -91.588;

n = -1.3;   Ke = (3.567*10^-12)*exp(90.13/8.314*T);
r = k*(P_A*(P_B)^2-(P_P/Ke))*(P_A*(P_B)^0.5)^n;

%differential equations
dF_AdW=-r;
dF_BdW=-2*r;
dF_PdW= r;
dTdW= (r*(-deltaH))/(F_A*CpA+F_B*CpB+F_V*CpV);

%Assemble Differential Equations
dUdW=[-r
    -2*r
    r
    (r*(-deltaH))/(F_A*CpA+F_B*CpB+F_V*CpV)];
end