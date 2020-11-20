% Physikalische Geodasesie Uebung 6
% Ziqing Yu 3218051
% 07/01/2020

clc
clear all
close all

rho = pi/180;
%% Task 1
l_num = 5;
theta = 0:180;
lambda = -180:180;
t_num = cosd(theta);
syms t l;
% Rodrigues-Ferrers
P_rf = zeros(length(theta),l_num + 1);
Y_rf = cell(l_num+1,2);

for i = 1:(l_num + 1)*2
    Y_rf(i) = {zeros(length(theta),length(lambda))};
end
P_rf_l = (t^2 - 1)^l;
for i = 1 : l_num
    P_rf_l = diff(P_rf_l,t);
end
P_rf_l = P_rf_l / (2^l_num * factorial(l_num));

% wenn m = 0, zonal spherical harmonics
% P
m = 0;
P_rf_lm = subs(P_rf_l, t, t_num);
P_rf_lm = double(subs(P_rf_lm, l, l_num));
P_rf_lm = P_rf_lm .* (1 - t_num.^2).^(m / 2);   
P_rf_lm = sqrt(2 * l_num + 1) * P_rf_lm;     
P_rf(:,1) = P_rf_lm'; % in 1. Spalte gespeichert.


% wenn m > 0
for m = 1 : l_num
    P_rf_l = diff(P_rf_l,t);
    P_rf_lm = subs(P_rf_l, l, l_num);
    P_rf_lm = subs(P_rf_lm, t, t_num);
    P_rf_lm = double(P_rf_lm);
    P_rf_lm = P_rf_lm .* (1 - t_num.^2).^(m / 2); 
    P_rf_lm = P_rf_lm * sqrt(2 * (2 * l_num + 1) * factorial(l_num - m) / factorial(l_num + m));
    P_rf(:,m + 1) = P_rf_lm';
end
figure
for i = 1 : l_num + 1 
    subplot(3,4,i)
    plot(90 - theta, P_rf(:,i))
    title(['m = ',num2str(i-1)])
    xlabel('\phi')
end
sgtitle('Rodrigues Ferres')


% Y
Y_1 = zeros(100,100);
Y_2 = zeros(100,100);
for m = 1:l_num+1
    Y_1 = P_rf(:,m) * cos(lambda * rho * (m-1));
    Y_2 = P_rf(:,m) * sin(lambda * rho * (m-1));
    Y_rf(m,1) = {Y_1};
    Y_rf(m,2) = {Y_2};
end

figure
for i = 1:l_num + 1
    subplot(3,4,i)
    hm = imagesc(cell2mat(Y_rf(i,1)));
    colormap('jet')
    xticks(0:90:360);
    xticklabels({'0','90','180','270','360'})
    yticks(-180:45:180);
    yticklabels({'90','45','0','-45','-90'})
    xlabel('\lambda, (-180 to 180)')
    ylabel('\theta, (0 to 180)')
    title(['m = ',num2str(i-1)])
end 
sgtitle('Using Rodriger Ferres')


% Recursive
P_rec = Normalized_Lengendre(l_num,theta);
figure
% jetzt das Ergebnis
for i = 1 : l_num + 1 
    subplot(3,4,i)
    plot(90 - theta, cell2mat(P_rec(l_num + 1,i)))
    title(['m = ',num2str(i-1)])
    xlabel('\phi')
end
sgtitle("Using recursive formulas")

% Y
Y_rec = cell(l_num+1,2);
for m = 1:l_num+1
    Y_1 = cell2mat(P_rec(l_num + 1,m))' * cos(lambda * rho * (m-1));
    Y_2 = cell2mat(P_rec(l_num + 1,m))' * sin(lambda * rho * (m-1));
    Y_rec(m,1) = {Y_1};
    Y_rec(m,2) = {Y_2};
end

% nur cos(lambda * m) soll gerechnet werden.
figure
for i = 1:l_num + 1
    subplot(3,4,i)
    hm = imagesc(cell2mat(Y_rec(i,1)));
    colormap('jet')
     xticks(0:90:360);
    xticklabels({'0','90','180','270','360'})
    yticks(-180:45:180);
    yticklabels({'90','45','0','-45','-90'})
    xlabel('\lambda, (-180 to 180)')
    ylabel('\theta, (0 to 180)')
    title(['m = ',num2str(i-1)])
end 
sgtitle('Using recursive formulas')
% xlabel und ylabel muss man noch einstellen!

%% Task 2
% wenn theta = 90, l = 0:100
l_task2 = 100;
theta90 = 90;
phi90 = 90 - theta90;
P_90 = Normalized_Lengendre(l_task2,theta90);

% wenn theta = 0:90, Loesungswerg ist aehnlich, aber komplizierter
theta090 = 0:90;
phi090 = 90 - theta090;
s090 = length(theta090);
Q_090 = Normalized_Lengendre(l_task2,theta090);

% jetzt rechnen wir die rechte Seite. 
task2_recht = cell(l_task2 + 1,1);
for il = 0 : l_task2
    % fuer jeder l gibt es ein 1 * 91 Matrix
    task2_recht(il + 1) = {zeros(1,s090)};
    for im = 0 : il
        task2_recht(il + 1) = {cell2mat(Q_090(il + 1,im + 1)) .* cell2mat(P_90(il + 1,im + 1)) + cell2mat(task2_recht(il + 1))};
    end
    task2_recht(il + 1) = {cell2mat(task2_recht(il + 1)) / (2 * il + 1)};
end

% linke Seite muss man, emmmmm wieder mit sym ableiten, weil da Legendre 
% Funktionen nicht nominiziert sind. 
psi = acosd(sind(phi90) .* sind(phi090) + cosd(phi90) .* cosd(phi090) .* cosd(zeros(1,length(phi90))));
t_num2 = cosd(psi);
l_num2 = 100;
task2_link = cell(l_num2 + 1,1);
P2 = (t.^2 - 1).^l;
P2_zwischen = subs(P2, t, t_num2);
P2_zwischen = double(subs(P2_zwischen,l,0));
task2_link(1) = {P2_zwischen * 1 / 1 / 1};
for i = 1 : l_num2
    P2 = diff(P2,t);
    P2_zwischen = subs(P2, t, t_num2);
    P2_zwischen = double(subs(P2_zwischen,l,i));
    task2_link(i+1) = {P2_zwischen / (2^i) / factorial(i)};
end


% Jetzt duerfen wir vergleichen.
diff_task2 = (cell2mat(task2_link) - cell2mat(task2_recht));
k = 9;
% und plotten
for i = 1 :fix((l_num2 + 1)/k)
    figure,hold on
    for j = 1 : k
        subplot(sqrt(k),sqrt(k),j)
        plot(diff_task2((i-1) * 9 + j,:))
    end
end

figure,hold on
for i = k*fix((l_num2 + 1)/k)+1:l_num2+1
    subplot(sqrt(k),sqrt(k),i-k*fix((l_num2 + 1)/k))
    plot(diff_task2(i,:))
end

%% Task3
l_task3 = 100;
theta3 = 0:180;
t3 = cosd(theta3);
P3 = Normalized_Lengendre(l_task3,theta3);

task3_recht = cell(l_task3 + 1,1);
for il = 0 : l_task3
    % fuer jeder l gibt es ein 1 * 181 Matrix
    task3_recht(il + 1) = {zeros(1,length(t3))};
    for im = 0 : il
        task3_recht(il + 1) = {cell2mat(P3(il + 1,im + 1)).^2 + cell2mat(task3_recht(il + 1))};
    end
    task3_recht(il + 1) = {cell2mat(task3_recht(il + 1)) / (2 * il + 1)};
end

for i = 1 :fix((l_task3 + 1)/k)
    figure,hold on
    for j = 1 : k
        subplot(sqrt(k),sqrt(k),j)
        plot(task3_recht{(i-1) * 9 + j})
    end
end

figure,hold on
for i = k*fix((l_task3 + 1)/k)+1:l_task3+1
    subplot(sqrt(k),sqrt(k),i-k*fix((l_task3 + 1)/k))
    plot(task3_recht{i})
end

%% task 4
fname = 'EGM96.txt';
fileid = fopen(fname,'r');

omega = 7.292115e-5;
% Datei einlesen
Cell = {};
i=1;
while ~feof(fileid)
    str=fgetl(fileid);
    if ~isempty(str)
        Cell{i} = str;
        a = split(str);
        degree(i, :) = str2double(a(2));
        order(i, :) = str2double(a(3));
        c_lm(i, :) = str2double(a(4));
        s_lm(i, :) = str2double(a(5));
        i=i+1;
    end
end

l4 = max(degree);
theta4 = 42 + 1;
lambda4 = 10 + 1;
r = 6379245.458;
R = 6378136.300;  
GM = 3.986004415e14;
t4 = cosd(theta4);
P4 = Normalized_Lengendre(l4, theta4);

Yc_4 = zeros(l4 + 1);
Ys_4 = zeros(l4 + 1);
clm = zeros(l4 + 1);
slm = zeros(l4 + 1);

for i=0:l4
    for j = 0:i
        clm(i+1,j+1) = c_lm(degree == i & order ==  j);
        Yc_4(i+1,j+1) = P4{i+1,j+1} * cosd((j) * lambda4);
        slm(i+1,j+1) = s_lm(degree == i & order ==  j);
        Ys_4(i+1,j+1) = P4{i+1,j+1} * sind((j) * lambda4);
    end
end
Vg = clm .* Yc_4 + slm .* Ys_4;
V_zwischen = sum(Vg,2);
Rdr = (R/r).^(1:l4+1);
V_gravitation = sum(V_zwischen .* Rdr') * GM / R;
V_centrifugal = 1 / 2 * omega^2 * r^2 * (sind(theta4))^2;
V_sum = V_gravitation + V_centrifugal;