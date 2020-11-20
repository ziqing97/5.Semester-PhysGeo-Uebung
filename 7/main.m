% Physikalischee Geodaesie Uebung 7
% Ziqing Yu 3218051
% 21/01/2020

clc
clear all;
close all,

Data = xlsread('Loop3.xlsx');
dh = Data(:,2);
ds = Data(:,3);
g = Data(:,4);
gama_45 = 9.806199203 * 10^5;
gama_mittel = 9.797644656 * 10^5;

%% task 1
sum_h = sum(dh);
g_mittel = (g(1:end-1) + g(2:end))/2;
sum_w = sum([0;g_mittel] .* dh(1:end));


%% task2
[H_d_45, H_O_45, DC_45, OC_45] = hoehe(gama_45,dh,g);
[H_d_mittel, H_O_mittel, DC_mittel, OC_mittel] = hoehe(gama_mittel,dh,g);

s = cumsum(ds);
figure, plot(s,H_d_45);
title('dynamic height')
figure, plot(s,H_O_45);
title('ortho height')
figure, plot(s, g);
title('gravity')
figure, plot(s,[0;DC_45]);
title('dynamic correction')
figure, plot(s,[0;OC_45]);
title('ortho correction')
figure, plot(s,[0;DC_mittel]);
title('dynamic correction with mean normal gravity')

function[H_d, H_O, DC, OC] = hoehe(gamma,dh,g)
%
H_O1 = 436.52;
H_d1 = (g(1) + 0.0424 * (H_O1))/gamma * H_O1;
%
dl = cumsum(dh);
dl = dl(2:end);
%
g_mittel = (g(1:end-1) + g(2:end))/2;
DC = cumsum((g_mittel - gamma)/gamma .* dh(2:end));
H_d = dl + DC + H_d1;
H_d = [H_d1;H_d];
%
H_O = gamma ./ (g(2:end) + 0.0424 * (H_O1 + dl)) .* H_d(2:end);
H_O = [H_O1;H_O];
%
OC = DC + (g(1) + 0.0424*H_O1 - gamma) / gamma * H_O1 - (g(2:end) + 0.0424 * H_O(2:end) - gamma) ./ gamma .* H_O(2:end);
end