function[a_x,a_y]=a_sphere(R,rho,X,Y,X_m,Y_m)
% thie function calculate the gravitational attraction of 
% a homogeneous solid sphere in 2D(z = 0).
% input: R: radius of the sphere
%        rho: density of the sphere
%        X,Y: coordinates of the points, where the attraction will be caculated, 
%             they can be vektor or matrix.
%        X_M,Y_M: coordinate of the center point of the sphere
% output: the gravitational attraction in x- and y- direction
    % constant 
    G = 6.672e-11;
    z = 0;
    % distance of the point and the centerpoint of the sphere
    r = sqrt((X - X_m).^2 + (Y - Y_m).^2 + z.^2);
    % 
    [x_size, y_size] = size(r);
    % calculate the value of the attraction
    a = zeros(x_size, y_size);
    a(r >= R) = -4 ./ 3 .* pi .* G .* rho .* R^3 ./ r(r >= R).^2;  
    a(r < R) = -4 ./ 3 .* pi .* G .* rho .* r(r < R); 
    a = -a; 
    % in vektor
    a_x = a .* (X - X_m) ./ r;
    a_y = a .* (Y - Y_m) ./ r;
end