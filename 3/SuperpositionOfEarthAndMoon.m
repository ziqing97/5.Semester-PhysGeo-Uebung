function [V, a_x, a_y] = SuperpositionOfEarthAndMoon(k)
% this function solve the task 3 and task 4, the superimposed gravitational
% field of the earth and moon will be visualized in the plane. 
    % constant
    R_E = 6371000;
    R_M = 1738000;
    m_E = 5.9736e24;
    m_M = 7.349e22;
    r_EM = 384400000;
    X_E = 0;
    Y_E = 0;
    X_M = r_EM * cosd(k * 10);
    Y_M = r_EM * sind(k * 10);
    rho_E = m_E / (4 / 3 * pi * R_E^3);
    rho_M = m_M / (4 / 3 * pi * R_M^3);
    % the range definition
    X = linspace(-r_EM / 2, 2 * r_EM, 100);
    Y = linspace(-r_EM / 2, 2 * r_EM, 100);   
    [x,y] = meshgrid(X, Y);
    % calculating the gravitational field
    V_E = V_sphere(R_E,rho_E,x,y,X_E,Y_E);
    V_M = V_sphere(R_M,rho_M,x,y,X_M,Y_M); 
    V = V_M + V_E;  
    V = V / 10^6;
    % presentation of the gravitational potential
    range = 0:0.1:2;
    figure, contour(X,Y,V,range,'ShowText','on');
    hold on 
    scatter(X_E, Y_E)
    scatter(X_M, Y_M)
    axis equal
    title('potential filed')
    % presentation of the gravitational attraction
    [aE_x,aE_y] = a_sphere(R_E,rho_E,x,y,X_E,Y_E);
    [aM_x,aM_y] = a_sphere(R_M,rho_M,x,y,X_M,Y_M);
    a_x = aE_x + aM_x;
    a_y = aE_y + aM_y;
    figure, quiver(x, y, a_x, a_y);
    hold on 
    scatter(X_E, Y_E)
    scatter(X_M, Y_M)
    title('vector filed')
end