function [V, a] = GravitationSphericalShell(r)
% this function solved the task 7, the gravitational potential and the
% attraction will be calculated and presented when the structure of the
% earth are core + mantle. 
% input r: the range (meter)
% output: V: the gravitational potential of r
%         a: the gravitational attraction of r
    % constant
    R_c = 3500000;
    R_m = 6400000;
    rho_c = 11200;
    rho_m = 4300;
    % calculate the potential
    V1 = V_shell(rho_c, r, R_c);
    V2 = V_shell(rho_m, r, R_c, R_m);
    V = V1 + V2;
    % present the potential
    figure, plot(r, V);
    title('potentil (0<r<4R)')
    xlabel('r, (m)')
    ylabel('potential (m^2/s^2)')
    % calculate the attraction
    a1 = a_shell(rho_c, r, R_c);
    a2 = a_shell(rho_m, r, R_c, R_m);
    a = a1 + a2;
    % present the attraction
    figure, plot(r, a);
    title('attraction (0<r<4R)')
    xlabel('r (m)')
    ylabel('attraction (g/s^2)')
end