% Physikalische Geodaesie Uebung 4 
% Ziqing Yu 3218051
% 03/12/2019
%% Initial
clc 
close all;
clear all;

%% constant
omega = 7.292115e-5; % rad/s
k = 1;

% Task 1 
% the latitude of the point.
phi_1 = (20 + 2 * k) / 180 *pi;
[V, V_c, W, a, g, dg, xi] = task1(phi_1);

% the disturbances within 0° < phi < 90°
phi_2 = linspace(0, pi/2);
[~, V_c2, ~, ~, ~, dg2, xi2] = task1(phi_2);

% plot everything we want
figure, plot(phi_2 / pi * 180, xi2 / pi * 180)
title('Disturbances of the direction')
xlabel('Latitude')
ylabel('\xi (degree)')

figure, plot(phi_2 / pi * 180, dg2)
title('Disturbances of the absolute value')
xlabel('Latitude')
ylabel('\sigma g (m/s^2)')

figure, plot(phi_2 / pi * 180, V_c2)
title('Centrifugal potential')
xlabel('Latitude')
ylabel('V_c (m^2/s^2)')



% Task 2
% constants
v = 400 / 3.6 / sqrt(2);    % to m/s und in NS und WE direction
phi_3 = 42 / 180 * pi;
theta = pi / 2 - phi_3;

% Coriolis acceleration
acor_ew = 2 * omega * [-cos(theta) * v; 0; sin(theta) * v];
acor_ns = 2 * omega * [0; cos(theta) * v; 0];

% Etoevoes correction
eto_ew = -acor_ew(3);
eto_ns = -acor_ns(3);

% Accuracy
eto_acc = 1e-5;
v_acc = eto_acc / (2 * omega * sin(theta));


