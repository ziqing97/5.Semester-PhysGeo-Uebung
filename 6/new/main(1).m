% Physikalische Geod?sie
% ?bung 6
%
% Junyang Gou
% 08.01.2020
%%
clear all
close all
clc
%% Initial values
l_max = 10;
theta = 0:1:180; % [deg]
t = cosd(theta);
t = t';

%% Legendre functions, recursive, the index is equal l+1 or m+1
% Initial values
Plm = cell(l_max+1);
Plm(:, :) = {zeros(size(t))};
% for i = 1:l_max+1
%     for j = 1:l_max+1
%         Plm{i, j} = zeros(size(t));
%     end
% end
Plm{0+1,0+1} = ones(size(t));
Plm{1+1,0+1} = t;
Plm{1+1,1+1} = sqrt(1 - t.^2);
% Recursive
for l = 2:l_max
    for m = 0:l
        if (m+1) == (l+1)
            Plm{l+1, m+1} = (2*l-1) * sqrt(1 - t.^2) .* Plm{(l-1)+1, (l-1)+1};
        else
            Plm{l+1, m+1} = 1/(l-m) * ((2*l-1) * t.*Plm{(l-1)+1, m+1} - ...
                (l-1+m) .* Plm{(l-2)+1, m+1});
        end
    end
end

%% Normalized Legendre functions, recursive,  the index is equal l+1 or m+1
% Initial values
Plm_norm = cell(l_max+1);
for i = 1:l_max+1
    for j = 1:l_max+1
        Plm_norm{i, j} = zeros(size(t));
    end
end
Plm_norm{0+1,0+1} = ones(size(t));
Plm_norm{1+1,0+1} = sqrt(3)*t;
Plm_norm{1+1,1+1} = sqrt(3*(1 - t.^2));
% Recursive
for l = 2:l_max
    for m = 0:l
        if (m+1) == (l+1)
            Plm_norm{l+1, m+1} = sqrt((2*l+1)/(2*l)) .* sqrt(1-t.^2) .*...
                Plm_norm{(l-1)+1, (m-1)+1};
        else
            Plm_norm{l+1, m+1} = sqrt((2*l+1)/((l+m)*(l-m))) * ...
                (sqrt(2*l-1) * t .* Plm_norm{(l-1)+1, m+1} - ...
                sqrt(((l-1+m)*(l-1-m))/(2*l-3)) .* Plm_norm{(l-2+1), m+1});
        end
    end
end

%% Plot the fully normalized legendre functions
figure
for i = 1:7
    plot(theta-90, Plm_norm{11, i});
    hold on
end
plot(theta-90, Plm_norm{11, 8}, 'r');
hold on
plot(theta-90, Plm_norm{11, 9}, 'g');
hold on
plot(theta-90, Plm_norm{11, 10}, 'b');
hold on
plot(theta-90, Plm_norm{11, 11}, 'y');
hold on

xlim([-90 90]);
legend
grid on

%% show the fully normalized spherical harmonics Ylm_norm
lambda = 0:1:360;
Ylm_norm_10 = cell(11, 1);
figure
for m = 0:10
    Ylm_norm_10{m+1, 1} = Plm_norm{11, m+1} * cosd(m*lambda);
    
    subplot(3,4,m+1)
    imagesc(Ylm_norm_10{m+1, 1})
    title_str = strcat('Y_{10,', num2str(m+1), '}(\theta, \lambda)');
    title(title_str);
    xlim([-10, 370]);
    ylim([-10, 190])
    xticks(0:40:360);
    yticks(0:20:180);
    colormap(jet);
    colorbar;
end

%% Task 2
theta_P = 90;
theta_Q = [0:1:90]';
phi_P = 90 - theta_P;
phi_Q = 90 - theta_Q';

cosPsi = sind(phi_P)*sind(phi_Q) + cosd(phi_P)*cosd(phi_Q);

syms l t
Pl_sym = 1/(2^l*factorial(l))*(t^2 - 1)^l;
Pl = cell(101, 1);
for i = 0:100
    diff_sym = diff(Pl_sym, t, i);
    Pl{i+1, 1} = double(subs(diff_sym, {t, l}, {cosPsi, i}))';
    disp(i);
end 

Plm_norm_P = Plm_recursive(100, theta_P, 'norm');
Plm_norm_Q = Plm_recursive(100, theta_Q, 'norm');

% Initial values
Plm_PQ = cell(101, 101);
right = cell(101, 1);
differenz = cell(101, 1);
for i = 1:101
    for j = 1:101
        Plm_PQ{i,j} = Plm_norm_P{i,j} .* Plm_norm_Q{i,j};
    end
    right{i, 1} = 1/(2*(i-1)+1) * sum(cell2mat(Plm_PQ(i,1:i)), 2);
    differenz{i, 1} = Pl{i, 1} - right{i, 1};
    disp(i);
end

% Plot the results
for i = 1:100
    if ~mod(i-1, 10)
        figure
    end
    subplot(5, 2, i-floor((i-1)/10)*10)
    plot(differenz{i, 1})
    title(i)
end

figure
plot(differenz{101, 1})
title(101)

%% Task 3
theta_3 = [0:180]';
Plm_3 = Plm_recursive(100, theta_3, 'norm');
right_3 = cell(101, 1);
for i = 1:101
    right_3{i, 1} = 1/(2*(i-1)+1) * sum(cell2mat(Plm_3(i,1:i)).^2, 2);
    disp(i);
end

%% Task 4
% Initial values
k = 6;
lambda_4 = 10+6;
theta_4 = 42+k;
r = 6379245.458;
R = 6378136.3;
GM = 3.986004415e14;
omega = 7.292115e-5;
% Import data
EGM96 = readEGM96('EGM96.txt');
clm = EGM96(:, 3);
slm = EGM96(:, 4);
%
Plm_4 = Plm_recursive(36, theta_4, 'norm');
V_temp_l = 0;
for l = 0:36
    V_temp_m = 0;
    for m = 0:l
        idx = EGM96(:, 1) == l & EGM96(:, 2) == m;
        c = clm(idx);
        s = slm(idx);
        V_temp_m = V_temp_m + Plm_4{l+1, m+1}*(c*cosd(m*lambda_4) + s*sind(m*lambda_4))
    end
    V_temp_l = V_temp_l + (R/r)^(l+1) * V_temp_m;
end

V = GM/R * V_temp_l;
Vc = 1/2*omega^2*r^2*(sind(theta_4))^2;
W = V + Vc;



