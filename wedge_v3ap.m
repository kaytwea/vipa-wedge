%   v3ap: air-spaced VIPA, point-sized waist, release version
%   - forumulated for constant wedge angle (alpha)
%   - variable image plane distance (zF)
%   - virtual source representation

%   Input parameters (can be hard-coded into function as desired)
%   n           [ ] index of refraction of air
%   R1          [ ] reflectivity of high reflector (intensity ratio)
%   R2          [ ] reflectivity of partial reflector
%   t0          [m] initial plate spacing of air-spaced VIPA
%   tilt        [°] tilt angle of VIPA
%   theta1/2    [°] min/max value of ray incident angles; based on input beam geometry
%   xF          [m] coordinate in imaging plane; variable 'v' in paper; must be row vector
%   numray      [ ] number of rays
%   pmax        [ ] number of internal reflections per ray
%   alpha       ["] wedge angle in arcseconds
%   d0          [m] distance along optical axis between VIPA pr and imaging lens
%   F           [m] focal length of imaging lens
%   zF          [m] image plane distance from lens; variable 'zi' in paper
%   k           [1/m] source wavenumber


function [I] = wedge_v3ap(n,R1,R2,t0,tilt,theta1,theta2,xF,numray,pmax,alpha,d0,F,zF,k)

r1 = sqrt(R1);          % reflection coeff of high reflector (reflected/incident amplitude ratio) 
r2 = sqrt(R2);          % reflection coeff of partial reflector
r1r2 = r1*r2;

alphad = alpha/3600;                        % [°] wedge angle in degrees

theta = linspace(theta1,theta2,numray);     % [°] ray incident angles

% ray tracking algorithm
tq = @(q,theta) (1-tand(theta-2*(q-1)*alphad)*tand(alphad))/(1+tand(theta-2*(q-1)*alphad)*tand(alphad));
rhoq = @(q,theta) tand(theta-2*(q-1)*alphad)+tand(theta-2*q*alphad);
phiq = @(q,theta) 1/cosd(theta-2*(q-1)*alphad)+1/cosd(theta-2*q*alphad);

I = zeros(length(k),length(xF));            % initializing fringe intensity pattern
R = r1r2.^(0:pmax);                         % row vector of field decay due to partial transmissions at second plate
R = repmat(R',[length(theta),length(xF)]);  % constructed for matrix multiplication

for g = 1:length(k)
    
    % del = path length difference [m] between source 0' and source p', at location xF along imaging plane; dominated by phi term (in original model)
    del = zeros(pmax+1,length(xF),length(theta));   % [source #, xF, incident angle]; row 1 is source 0
    del(1,:,:) = repmat(-1/(2*zF)*d0*xF.^2/(zF+d0*(1-zF/F)),[1,1,length(theta)]);    % initializing source 0 to phase term when rho = phi = 0
    
    for j = 1:length(theta)         % cycle through rays

                                    % t, rho, phi all row vectors for source 1 to pmax; column 1 is source 1
        t = NaN(1,pmax);            % plate thickness [m] at source p
        rho = NaN(1,pmax);          % distance [m] along second mirror surface from source 0 to p
        phi = NaN(1,pmax);          % path length [m] of ray from source 0 to p
        
        voffset = NaN(1,pmax);      % [m] vertical offset of virtual source
        zoffset = NaN(1,pmax);      % [m] horizontal offset of virtual source
        phiv = NaN(1,pmax);         % phase offset of virtual source, i.e. 0
        
        t_seg = zeros(1,pmax);      % segment t reduction scalar; column 1 is source 1
        rho_seg = zeros(1,pmax);    % segment rho; col 1 is source 1
        phi_seg = zeros(1,pmax);    % segment phi; col 1 is source 1

        for i = 1:pmax              % cycle through internal reflections
                                    % compile plate thickness, rho, phi from source 1 to pmax; then convert to virtual source locations

            t_seg(i) = tq(i,theta(j));
            rho_seg(i) = rhoq(i,theta(j));
            phi_seg(i) = phiq(i,theta(j));
            
            t(i) = t0*prod(t_seg(1:i));
            rho(i) = sum(t(1:i).*rho_seg(1:i));
            phi(i) = sum(t(1:i).*phi_seg(1:i));

            voffset(i) = +rho(i)*cosd(tilt) + phi(i)*sind(tilt-theta(j)+2*i*alphad);
            zoffset(i) = -rho(i)*sind(tilt) + phi(i)*cosd(tilt-theta(j)+2*i*alphad);
            phiv(i) = 0;
            
        end

        del(2:end,:,j) = repmat(+zoffset',[1,length(xF)]) ...
                        + 1/2*(1-zF/F)*repmat(voffset'.^2,[1,length(xF)])./repmat(zF+(d0+zoffset')*(1-zF/F),[1,length(xF)]) ...
                        - voffset'*xF./repmat(zF+(d0+zoffset')*(1-zF/F),[1,length(xF)]) ...
                        - 1/(2*zF)*(d0+zoffset')*xF.^2./repmat(zF+(d0+zoffset')*(1-zF/F),[1,length(xF)]) ...
                        + repmat(phiv',[1,length(xF)]);

    end % theta cycle

    del = permute(del,[1 3 2]);                                 % permute to [source #, incidence angle, xF] for reshaping
    del = reshape(del,[size(del,1)*size(del,2),size(del,3)]);   % reshapes to 2D matrix; each column is all sources for every incident angle at single xF
                                                                % eg. column 1 is source 0 to pmax for angle 1, source 0 to pmax for angle 2, etc, at xF(1)
    
    % computation of intensity profile across xF for wavenumber k; field amplitude summation across sources and incident angle
    I(g,:) = (sum(R.*cos(2*pi*k(g)*n*del),'omitnan')).^2 + (sum(R.*sin(2*pi*k(g)*n*del),'omitnan')).^2;

end % k cycle

end % function