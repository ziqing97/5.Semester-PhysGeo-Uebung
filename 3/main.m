% Phisikalische Geodaesie Uebung 3
% Ziqing Yu 3218051
% erstellt am 19/11/2019
clc
clear all;
close all;
% some constants to be needed.
G = 6.672e-11;
k = 1;

%% Superposition of Earth's and Moon's gravitational fields
% task 2 to 4
[V, a_x, a_y] = SuperpositionOfEarthAndMoon(k);

%% Gravitational potential and attraction of spherical shells
% task 7 
% constant
R_m = 6400000;
% definite the range
r_2 = linspace(0, 4 * R_m);
% calculate and present
[V_earth, a_earth] = GravitationSphericalShell(r_2);

%% PREM density model of the Earth
% task 8
% read the data
load PREM.mat
% definite the range 
r_3 = linspace(0, 2 * 6371000);
% calculate the potential and the attraction
[V_PREM, a_PREM] = PREM_densityModel(PREM, r_3);
% presentation
figure, plot(r_3, V_PREM)
title('gravitational potential of the PREM model')
xlabel('r, (m)')
ylabel('potential (m^2/s^2)')
figure, plot(r_3, a_PREM)
title('gravitational attraction of the PREM model')
xlabel('r (m)')
ylabel('attraction (g/s^2)')
% task 9 ,the potential and the attraction of the position Earth surface
R_E = 6371000;
[V_sur, a_sur] = PREM_densityModel(PREM, R_E);
% task 10
[~, posmax_V] = max(V_PREM);
[~, posmax_a] = max(a_PREM);
% the position and the value of the max. potential
rmax_V = r_3(posmax_V);
V_max = V_PREM(posmax_V);
% the position and the value of the max. attraction
rmax_a = r_3(posmax_a);
a_max = a_PREM(posmax_a);
