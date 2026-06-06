% Code Name: ASG2 Variable OPR BPR
% Code Description: ASG2 Variable OPR and BPR trend calculator and plotter for turbofan engine
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE403 - Section 03DB
% Date: 02/25/2026
% Worked With: Copilot.
%Iniciation
clc; clear; close all;

% Givens
M0 = 0.8;
T0 = 218.8; 
P0 = 23.8e3; 
R  = 287;
gamma_c = 1.4;
cp_c    = 1004;
gamma_h = 1.33;
cp_h    = 1152;
FPR = 1.5;
Df  = 3;
M2  = 0.5;
TET = 1700;
eta_c = 0.88;
eta_t = 0.90;
eta_m = 0.99;
pi_d  = 0.98;
pi_b  = 0.95;
pi_bp = 0.97;
Cv = 0.99;
LHV = 43.15e6;
g0  = 9.81;

% Freestream and fan calculations
a0  = sqrt(gamma_c*R*T0);
V0  = M0*a0;
Tt0 = T0*(1 + (gamma_c-1)/2*M0^2);
Pt0 = P0*(1 + (gamma_c-1)/2*M0^2)^(gamma_c/(gamma_c-1));

Pt2 = Pt0*pi_d;
Tt2 = Tt0;
T2 = Tt2/(1+(gamma_c-1)/2*M2^2);
P2 = Pt2/(1+(gamma_c-1)/2*M2^2)^(gamma_c/(gamma_c-1));
rho2 = P2/(R*T2);
A_fan = pi*(Df/2)^2;
V2 = M2*sqrt(gamma_c*R*T2);
mdot_total = rho2*A_fan*V2;

% Thermal Efficiency versus varibale OPR
BPR_fixed = 10;
OPR_vec = [30,40,50,60,70,80];

mdot_core = mdot_total/(1+BPR_fixed);
mdot_bp   = mdot_total - mdot_core;

PR_HPC = OPR_vec ./ FPR;
Tt3 = Tt2 .* (1 + (PR_HPC.^((gamma_c-1)/gamma_c)-1)/eta_c);
Pt3 = Pt2 .* PR_HPC;

Tt4 = TET;
f = (cp_h*Tt4 - cp_c.*Tt3) ./ (LHV - cp_h*Tt4);
mdot_f = mdot_core .* f;
W_HPC = mdot_core .* cp_c .* (Tt3 - Tt2)*eta_t;

Tt13 = Tt2*(1 + (FPR^((gamma_c-1)/gamma_c)-1)/eta_c);
W_fan = mdot_total*cp_c*(Tt13 - Tt2)*eta_t;

Tt45 = Tt4 - W_HPC./(mdot_core*cp_h*eta_m);
Tt5  = Tt45 - W_fan./(mdot_core*cp_h*eta_m);

Ve_core = sqrt(2*Cv*cp_h.*Tt5);
Ve_bp   = sqrt(2*Cv*cp_c*Tt13);

JetPower = 0.5*mdot_core.*(1+f).*Ve_core.^2 + 0.5*mdot_bp*Ve_bp^2 - 0.5*mdot_total*V0^2;
thermal_eff = ((JetPower/2) ./ ((mdot_f)*LHV))-1;

% Proposive Efficiency with variable BPR
OPR_fixed = 50;
BPR_vec = [2,5,10,15];

mdot_core = mdot_total ./ (1+BPR_vec);
mdot_bp   = mdot_total - mdot_core;

PR_HPC = OPR_fixed/FPR;
Tt3 = Tt2*(1 + (PR_HPC^((gamma_c-1)/gamma_c)-1)/eta_c);

f = (cp_h*Tt4 - cp_c*Tt3) / (LHV - cp_h*Tt4);
mdot_f = mdot_core .* f;

W_HPC = mdot_core .* cp_c .* (Tt3 - Tt2);
W_fan = mdot_total*cp_c*(Tt13 - Tt2);

Tt45 = Tt4 - W_HPC./(mdot_core*cp_h*eta_m);
Tt5  = Tt45 - W_fan./(mdot_core*cp_h*eta_m);

Ve_core = sqrt(2*Cv*cp_h.*Tt5);
Ve_bp   = sqrt(2*Cv*cp_c*Tt13);

F_core = mdot_core.*(1+f).*Ve_core;
F_bp   = mdot_bp.*Ve_bp;
F_net  = F_core + F_bp - mdot_total*V0;

JetPower = 0.5*mdot_core.*(1+f).*Ve_core.^2 + 0.5*mdot_bp.*Ve_bp.^2 - 0.5*mdot_total*V0^2;
prop_eff = F_net*V0 ./ JetPower;

% Plots
figure;
plot(OPR_vec, thermal_eff,'LineWidth',2)
grid on
xlabel('Overall Pressure Ratio (OPR)')
ylabel('Thermal Efficiency')
title('Thermal Efficiency vs OPR')

figure;
plot(BPR_vec, prop_eff,'LineWidth',2)
grid on
xlabel('Bypass Ratio (BPR)')
ylabel('Propulsive Efficiency')
title('Propulsive Efficiency vs BPR')
