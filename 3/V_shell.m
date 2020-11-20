function[V] = V_shell(rho, r, R_in, R_out)
% this function calculate the gravitational potential of a homogeneous
% shell or sphere
% input: rho: the density of the shell
%        r: the distance of the shellcenter and the point, where the
%           potential will be calculated. it can be a vector or a matrix
%        R_in: the radius of the inner radius or the radius of the sphere
%        R_out: the radius of the outer radius
% output: V: the gravitational potential

    % constant
    G = 6.672e-11;
    % if the nargin is 4, then the potential of a shell will be
    % calculated, if the nargin is 3, then the potential of a sphere will
    % be calculated. any other cases will be regarded as error. 
    if nargin == 4
        V = 2 * pi * G * rho * (R_out^2 - r.^2/3) - 4 / 3 * pi * G * rho * R_in^3 ./ r;
        V(r < R_in) = 2 * pi * G * rho * (R_out^2 - R_in^2);
        V(r >= R_out) = 4 / 3 * pi * G * rho *(R_out^3 - R_in^3) ./ r(r >= R_out);
    elseif nargin == 3
        R = R_in;
        V(r < R) = 2 * pi * G * rho * (R^2 - r(r < R).^2 ./ 3);
        V(r > R) = 4 / 3 * pi * G * rho * R^3 ./ r(r > R);
    else
        error('wrong with input')
    end
end