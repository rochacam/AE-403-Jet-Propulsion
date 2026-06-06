% Code Name: Nozzle Inlet Mach Number Interator
% Code Description: Iterates the mass flow function to find the mach number for the inlet of a
% conv-div nozzle
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE403 - Section 03DB
% Date: 01/26/2026
% Worked With: N/a

% Iniciation
clear; clc;

% Given data and constant
gamma_after = 1.33;
P07 = 0.99 * 850952.8;
Pa = 30090;
A7 = 0.0330844;
R = 287;
T07 = 2700 - 50; 
mass_flow = 13.5449;

% Setting Up Newtonian Interation
Const = (P07*A7)/sqrt(T07)*sqrt(gamma_after/R);
target_m = mass_flow;
f_m = @(M) Const*M*(1+(gamma_after-1)/2*M^2)^(-(gamma_after+1)/(2*(gamma_after-1)))-target_m;

% Initial guess for Mach number
M_guess = 0.1;
M7 = fzero(f_m,M_guess);

% Output Results
fprintf('Calculated Mach Number (M7) : %.4f\n' , M7);

