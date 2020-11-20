function[V,r]=V_sphere(R,rho,X,Y,X_m,Y_m)
% thie function calculate the gravitational potential of 
% a homogeneous solid sphere in 2D (z = 0).
% input: R: radius of the sphere
%        rho: density of the sphere
%        X,Y: coordinates of the points, where the potential will be caculated, 
%             they can be Vektor or Matrix.
%        X_M,Y_M: coordinate of the center point of the sphere
% output: the gravitational potential
    % constant
    G = 6.672e-11;
    z = 0;
    % distance of the point and the centerpoint of the sphere
    r = sqrt((X - X_m).^2 + (Y - Y_m).^2 + z.^2);
    [x_size, y_size] = size(r);
    % caculate the potential
    V = zeros(x_size, y_size);
    V(r >= R) = 4 ./ 3 .* pi .* G .* rho .* R^3 ./ r(r >= R);
    V(r < R) = 2 .* pi .* G .* rho .* (R.^2 - (r(r < R).^2) ./ 3);
end