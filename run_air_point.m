%   This script runs the air-spaced VIPA model with input beam converging to point-sized waist: wedge_v3ap

n = 1;          % index of refraction of air

R1 = 0.995;     % reflectivity of high reflector (intensity ratio)
R2 = 0.950;     % reflectivity of partial reflector

d0 = 0.10;      % [m] approx distance between VIPA and imaging lens; output relatively insensitive to this
t0 = 15.0e-3;   % [m] initial plate spacing of air-spaced VIPA

%   Computational parameters - these have significant effect on computation time
%   Note: for air-spaced + point waist ONLY, numray = 1 sufficient as all virtual sources are co-located
%   -----------------------------------------------------------------------
k = linspace(18784.90e2,18785.20e2,13);     % [1/m] wavenumber of fringe patterns to compute
xF = (5:.001:30)*1e-3;                      % [m] range and step size of computational domain (i.e. image plane)
numray = 1;                                 % number of rays
pmax = 200;                                 % number of internal reflections per ray
%   -----------------------------------------------------------------------

F = 2000e-3;        % [m] focal length of imaging lens
alpha = 0.00;       % ["] wedge angle in arcseconds
zF = 2000e-3;       % [m] image plane distance from lens

% Fc = 300e-3;            % [m] focal length of cylindrical lens
% beam_diam = 5e-3;       % [m] diameter of collimated beam incident on cylindrical lens
tilt = 0.10;        % [°] tilt angle of VIPA
theta1 = 0.000;     % [°] lower limit of incident angle range
theta2 = 0.577;     % [°] upper limit of incident angle range

tic
[I] = wedge_v3ap(n,R1,R2,t0,tilt,theta1,theta2,xF,numray,pmax,alpha,d0,F,zF,k);
toc

figure()
plot(xF*1e3,I)
title({"air-spaced VIPA, input beam w/ point-sized waist" ...
    "\theta = " + theta1 + "° - " + theta2 + "°, tilt = " + tilt + "°, \alpha = " + alpha + " arcsec, zF = " + zF*1e3 + " mm" ...
    "numray = " + numray + ", p = " + pmax + ", domain resolution: " + (xF(2)-xF(1))*1e6 + " µm"})
legend(string(k'/100))
legend off
xlabel('Image plane [mm]')
ylabel('Fringe intensity [a.u.]')



