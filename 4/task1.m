function[V, V_c, W, a_abs, g_abs, dg, xi] = task1(phi)
% This Function solve the first task, the input phi is the latitude, which
% can be a vector. the output V is the gravitational potential, V_c is the
% centrifugal potential, W is the gravity potential, a_abs is the absolute
% value of teh gravitational attraction, g_abs is the absolute value of the
% gravity attraction, dg is the disturbance of the attraction g-a and xi is
% the disturbance of the direction. All of the outputs can be vectors with
% the same size of phi

% constants
rho = 5515;      % kg/m^3
omega = 7.292115e-5; % rad/s
lambda = 10 / 180 * pi;
r = 6371000;
R = 6371000;
G = 6.672e-11;

% potential
V = 4 / 3 * pi * G * rho * R^3 / r;
V_c = 1 / 2 *omega^2 * r^2 * (cos(phi)).^2;
W = V + V_c;

% attraction
x = r * cos(phi) * cos(lambda);
y = r * cos(phi) * sin(lambda);
z = r * sin(phi);
len = length(x);
a = -4 / 3 / r^3 * pi * G * rho * R^3 * [x; y; z];
a_c = omega^2 * [x; y,; zeros(1,len)];
g = a + a_c;

% absolute value
a_abs = sqrt(sum(a.^2,1));
g_abs = sqrt(sum(g.^2,1));

% the direction and the disturbance
parameter1 = (a .* g);
parameter2 = sum(parameter1,1);
xi = acos( parameter2 ./ (a_abs .* g_abs));
dg = g_abs - a_abs;
end