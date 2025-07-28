%   This script runs the solid-body VIPA model with input beam converging to point-sized waist: wedge_v3sp

n = 1.4;        % index of refraction of solid-body VIPA

R1 = 0.995;     % reflectivity of high reflector (intensity ratio)
R2 = 0.965;     % reflectivity of partial reflector

d0 = 0.10;      % [m] approx distance between VIPA and imaging lens; output relatively insensitive to this
t0 = 5.0e-3;    % [m] initial mirror spacing of solid-body VIPA

%   Computational parameters - these have significant effect on computation time
%   Note: due to computational load of 1000 rays and 200 internal reflections, domain is subdivided
%   -----------------------------------------------------------------------
k = linspace(18784.90e2,18785.20e2,13);     % [1/m] wavenumber of fringe patterns to compute
xF = (5:.005:30)*1e-3;                      % [m] range and step size of computational domain (i.e. image plane)
    xF(end) = [];                                % delete last entry for even number of coordinates
    subnum = 20;                                 % number of domain subdivisions
    xFsub = reshape(xF,length(xF)/subnum,[]);    % length(xF)/subnum = length of subdomain; must be integer
    xFsub = xFsub';
numray = 3000;                              % number of rays
pmax = 200;                                 % number of internal reflections per ray
%   -----------------------------------------------------------------------

F = 2000e-3;            % [m] focal length of imaging lens
alpha = 0.00;           % ["] wedge angle in arcseconds
zF = 2000e-3;           % [m] image plane distance from lens

% Fc = 300e-3;            % [m] focal length of cylindrical lens
% beam_diam = 5e-3;       % [m] diameter of collimated beam incident on cylindrical lens
tilt = 0.10;            % [°] tilt angle of VIPA
theta1 = 0.000;         % [°] lower limit of incident angle range
theta2 = 0.577+3;       % [°] upper limit of incident angle range
                        % output improved by artificial increase to theta range; 3° added to theta2

tic
I = [];                 % full fringe pattern
for i = 1:subnum        % crunch fringe pattern fragments in sequence
    [Isub] = wedge_v3sp(n,R1,R2,t0,tilt,theta1,theta2,xFsub(i,:),numray,pmax,alpha,d0,F,zF,k);
    I = cat(2,I,Isub);  % stitch fringe pattern segment to full pattern
end
toc

figure()
plot(xF*1e3,I)
title({"solid-body VIPA, n = " + n + ", input beam w/ point-sized waist" ...
    "\theta = " + theta1 + "° - " + theta2 + "°, tilt = " + tilt + "°, \alpha = " + alpha + " arcsec, zF = " + zF*1e3 + " mm" ...
    "numray = " + numray + ", p = " + pmax + ", domain resolution: " + (xF(2)-xF(1))*1e6 + " µm"})
legend(string(k'/100))
legend off
xlabel('Image plane [mm]')
ylabel('Fringe intensity [a.u.]')
