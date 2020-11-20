function [s_V, s_a] = PREM_densityModel(Model, r )
% this function solve the task 8 and task 9, the potential and the
% attraction will be calculated when the PREM density model of the earth
% will be used. the input Model are the PREM model(the first colume is the radius
% the second is the density), the input r is the range. the output s_V and
% s_a are the gravitation potential and attraction of the earth. 
    R = Model(:,1) * 1000; % km to m
    rho = Model(:,2);
    len1 = length(R);
    len2 = length(r);
    V = zeros(len1,len2);
    a = zeros(len1,len2);
    % calculate the P and a of the first layer
    V(1,:) = V_shell(rho(1), r, R(1));
    a(1,:) = a_shell(rho(1), r, R(1));
    % calculate the P and a of the secoud to the last layer
    for i = 2 : len1
        V(i,:) = V_shell(rho(i), r, R(i - 1), R(i));
        a(i,:) = a_shell(rho(i), r, R(i - 1), R(i));
    end
    % the total potential and the attraction
    s_V = sum(V, 1);
    s_a = sum(a, 1);
end