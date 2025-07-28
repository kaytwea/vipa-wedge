%   v3ab: air-spaced VIPA, finite beam waist, release version
%   - forumulated for constant wedge angle, alpha
%   - variable image plane distance, zF
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
%   beam_diam   [m] diameter of collimated beam incident upon cylindrical lens
%   Fc          [m] focal length of cylindrical lens
%   d0          [m] distance along optical axis between VIPA pr and imaging lens
%   F           [m] focal length of imaging lens
%   zF          [m] image plane distance from lens; variable 'zi' in paper
%   k           [1/m] source wavenumber


function [I] = wedge_v3ab(n,R1,R2,t0,tilt,theta1,theta2,xF,numray,pmax,alpha,beam_diam,Fc,d0,F,zF,k)

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

    w0 = (1/k(g))*Fc/(pi*beam_diam/2);                          % [m] radius of (Gaussian) beam waist 
    w0x = linspace((tilt-theta1)/(theta2-tilt)*w0,-w0,numray);  % [m] radius of specific ray at beam waist; pre-factor for w0 accounts for rays that contribute (in theory) to beam waist, but fail to couple into VIPA
                                                                % add offset to w0x to implement VIPA height adjustment; no accounting for any associated input clipping
                                                                    
    gauss = @(x) exp(-x.^2/w0^2);           % Gaussian field (not intensity) profile
    beam  = gauss(w0x);                     % beam field profile at waist
    beamop = repmat(beam,pmax+1,1);
    beamop = repmat(beamop(:),1,length(xF));
    
    % del = path length difference [m] between source 0' and source p', at location xF along imaging plane; dominated by phi term (in original model)
    del = zeros(pmax+1,length(xF),length(theta));   % [source #, xF, incident angle]; row 1 is source 0
    
    for j = 1:length(theta)         % cycle through rays at variable incidence theta

                                    % t, rho, phi all row vectors for source 0 to pmax; column 1 is source 0
        t = NaN(1,pmax+1);          % [m] plate thickness at source p
        rho = NaN(1,pmax+1);        % [m] distance along second mirror surface from source 0 to p (+ initial offset from beam waist geometry)
        phi = NaN(1,pmax+1);        % [m] path length of ray from source 0 to p (+ initial offset from beam waist geometry)
       
        voffset = NaN(1,pmax+1);    % [m] vertical offset of virtual source
        hoffset = NaN(1,pmax+1);    % [m] horizontal offset of virtual source
        phiv = NaN(1,pmax+1);       % phase offset of virtual source, i.e. phi(1)
        
        t_seg = zeros(1,pmax+1);    % segment t reduction scalar; col 1 is source 0
        rho_seg = zeros(1,pmax+1);  % segment rho; col 1 is source 0
        phi_seg = zeros(1,pmax+1);  % segment phi; col 1 is source 0
        
        t_seg(1) = 1;                       % NOTE: for coherent beam, coincidence of beam waist and mirror surface is considered origin - optical axis if no offset to w0x
        rho_seg(1) = w0x(j)/cosd(tilt);     % hypotenuse of radius at waist
        phi_seg(1) = w0x(j)*tand(tilt);     % additional z propagation

        t(1) = t0*t_seg(1);                 % source 0: assume constant t(1) across beam waist
        rho(1) = rho_seg(1);                % source 0: initial offset along partial reflector surface
        phi(1) = phi_seg(1);                % source 0: initial phase offset from additional z propagation
        
        voffset(1) = w0x(j);                % "virtual" source 0 is vertical beam waist
        hoffset(1) = 0;
        phiv(1)    = 0;
        
        for i = 2:pmax+1                    % compile plate thickness, rho, phi from source 1 to pmax
            
            t_seg(i) = tq(i-1,theta(j));
            rho_seg(i) = rhoq(i-1,theta(j));
            phi_seg(i) = phiq(i-1,theta(j));
            
            t(i) = t(1)*prod(t_seg(2:i));
            rho(i) = sum(t(2:i).*rho_seg(2:i)) + rho(1);
            phi(i) = sum(t(2:i).*phi_seg(2:i)) + phi(1);
                        
            voffset(i) = +rho(i)*cosd(tilt) + phi(i)*sind(tilt-theta(j)+2*(i-1)*alphad);
            hoffset(i) = -rho(i)*sind(tilt) + phi(i)*cosd(tilt-theta(j)+2*(i-1)*alphad);
            phiv(i) =    phiv(1);

        end
        
        del(1:end,:,j) = repmat(+hoffset',[1,length(xF)]) ...
                        + 1/2*(1-zF/F)*repmat(voffset'.^2,[1,length(xF)])./repmat(zF+(d0+hoffset')*(1-zF/F),[1,length(xF)]) ...
                        - voffset'*xF./repmat(zF+(d0+hoffset')*(1-zF/F),[1,length(xF)]) ...
                        - 1/(2*zF)*(d0+hoffset')*xF.^2./repmat(zF+(d0+hoffset')*(1-zF/F),[1,length(xF)]) ...
                        + repmat(phiv',[1,length(xF)]);

    end % theta cycle
        
    del = permute(del,[1 3 2]);                                 % permute to [source #, incidence angle, xF] for reshaping
    del = reshape(del,[size(del,1)*size(del,2),size(del,3)]);   % reshapes to 2D matrix; each column is all sources for every incident angle at single xF
                                                                % eg. column 1 is source 0 to pmax for angle 1, source 0 to pmax for angle 2, etc, at xF(1)
    
    % computation of intensity profile across xF for wavenumber k; field amplitude summation across sources and incident angle
    I(g,:) = (sum(R.*beamop.*cos(2*pi*k(g)*n*del),'omitnan')).^2 + (sum(R.*beamop.*sin(2*pi*k(g)*n*del),'omitnan')).^2;
%     I(g,:) = (sum(R.*cos(2*pi*k(g)*n*del),'omitnan')).^2 + (sum(R.*sin(2*pi*k(g)*n*del),'omitnan')).^2;     % no Gaussian beam profile

end % k cycle

end % function