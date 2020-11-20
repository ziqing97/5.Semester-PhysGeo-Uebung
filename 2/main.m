% Phisikalische Geadaesie Uebung 2
% Ziqing Yu
% 3218051
clc;
clear all;
close all;

% Datei einlesen
fname = 'B2100200GGP.sec';
fileid = fopen(fname,'r');

Cell = {};
i=1;
while ~feof(fileid)
    str=fgetl(fileid);
    if ~isempty(str)
        if str(1) == '2' 
            Cell{i} = str;
            spli = split(str);
            day(i, :) = (spli(1));
            hour(i, :) = (spli(2));
            time(i, :)=str2double(strcat(day(i,:),hour(i,:)));
            grav(i, :) = str2double(spli(3));
            pressure(i, :) = str2double(spli(4));
            i=i+1;
        end
    end
end


% fast foureir transformation
grav_ugal=grav*(-80.006);
plot(grav_ugal)
xticks([find(time==20100201000000),find(time==20100211000000),find(time==20100221000000), ...
    find(time==20100228235900)]);
xticklabels({'1','11','21','1 March'});
ylim([-500,50])
xlabel('February 2010')
ylabel('Gravity [\muGal]')



figure 
grav_f=fft(grav_ugal);
% timescalar=1 min, f=1/60 HZ
FS=1/60;
L=length(grav_f);
P2=abs(grav_f/L);
P1=P2((1:L/2+1));
P1(2:end-1)=2*P1(2:end-1);
f=FS*(0:(L/2))/L;
plot(f,P1);
ylabel('Gravity')
xlabel('Frequenz [HZ]')




