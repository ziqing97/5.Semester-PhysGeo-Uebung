% Physikalische Geodaesie Uebung 5
% Ziqing Yu
% 3218051

clc
clear all;
close all;

load loisach.mat;
G = 6.672e-11;
%% Task 1
dx = 100;
dy = 100;
gr = gr + 98;
dM1 = sum(sum(gr)) * 10^(-5);
dM1 = dM1 * dx * dy / (4 * pi * G);

%% Task 2
dM2 = 300 * 2000 * 5000 * (-700) / 2;












