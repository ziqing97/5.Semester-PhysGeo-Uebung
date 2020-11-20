clc
clear all;
close all;
%% gegeben Daten
h13=6.897; % m
h16=17.26; % m
h36=h16-h13;

GM=3.986004415e14; % m^3/s^2
R=6378136; % m

sigma_h=0.03; % m
fname = 'Messergebnis.txt';
fileid = fopen(fname,'r');

Cell = {};
i=1;
while ~feof(fileid)
    str=fgetl(fileid);
    if ~isempty(str)
        if str(2) =='1' || str(2)=='2'
            Cell{i} = str;
            
            a = split(str);
            grav(i, :) = str2double(a(5));
            dur(i, :) = str2double(a(11));
            i=i+1;
        end
    end
end

ob=zeros(6,1);
for i=1:6
    ob(i)=(grav(i*3)+grav(i*3-1)+grav(i*3-2))/3;
end
l=zeros(5,1);
for i=1:5
    l(i)=ob(i+1)-ob(i);
end

A1=[1;-1;0;1;-1];
A2=[0;0;1;-1;1];
A3=zeros(5,1);
A=[A1,A2,A3];
x=pinv(A'*A)*A'*l;



