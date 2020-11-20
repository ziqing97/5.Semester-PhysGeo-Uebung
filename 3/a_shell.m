function[a] = a_shell(rho, r, R_in, R_out)
% this function calculate the gravitational attraction of a homogeneous
% shell or sphere
% input: rho: the density of the shell
%        r: the distance of the shellcenter and the point, where the
%           attraction will be calculated. it can be a vector or a matrix
%        R_in: the radius of the inner radius or the radius of the sphere
%        R_out: the radius of the outer radius
% output: a: the gravitational attraction

    % constant
    G = 6.672e-11;
    % if the nargin is 4, then the attraction of a shell will be
    % calculated, if the nargin is 3, then the attraction of a sphere will
    % be calculated. any other cases will be regarded as error.
    if nargin == 4
        a = -4 / 3 * pi* G * rho * (r.^3 - R_in^3) ./ r.^2;
        a(r < R_in) = 0;
        a(r >= R_out) = -4 / 3 * pi * G * rho * (R_out^3 - R_in^3) ./ r(r >= R_out).^2;
    elseif nargin == 3
        R = R_in;
        a(r < R) = -4 / 3 * pi * G * rho .* r(r < R);
        a(r > R) = -4 / 3 * pi * G * rho * R^3 ./ r(r > R).^2;
    else 
        error('wrong with input')
    end
    a = -a;
end