function [ flp ] = p_art_water( flp, sph, geom )

% function [ flp ] = p_art_water( flp, sph, geom )
% Computes pressure and sound speed for artificial water.
% Created:     2011
% Last change: 24.06.2021

% %==========================================================================
% % Equation of state for artificial water
% % Form 1 (Cole, Batchelor 1973, Monaghan)
% %==========================================================================
% % Calculating B, limit for the maximum density variation [Pa].
% % If B has the following form:
% %    B = 100*rho_0*v^2/gamma,
% % then the relative density fluctuactions should be about 1% (i.e., the
% % water will be treated as a very slightly compressible fluid).
% % In the above expression, v is the maximum velocity. According to Monaghan
% % (1994), "when a dam of height H collapses, an appropriate upper bound to
% % the speed of water is given by 2*g*H." Monaghan and Kos (1999) tell to
% % use g*H instead. See alos page 65 of my first research notebook.
% 
% v_max = sqrt( -2 * sph.g * (geom.y_max - geom.y_min) );
% gamma = 7;
% rho_0 = 1000;        % reference density [kg/m^3]
% B = 100 * rho_0 * v_max^2/gamma;        
% % Atmosferic pressure = 101300 Pa
% % B = 101300;
% 
% % Equation of state for artificial water (Cole, 1948; Batchelor, 1974)
% flp.p = B * ( (flp.rho./rho_0).^gamma - 1 );
% 
% % Speed of sound [m/s] for artificial water
% % speed of sound in freshwater at 25°C: 1497 m/s
% % speed of sound in seawater: 1560 m/s
% flp.c = 10 * v_max * ones( geom.tnp, 1 );
% % c = 1480;

%==========================================================================
% Equation of state for artificial water
% Form 2 (Morris, 1997)
%==========================================================================
% Define the speed of sound for the total number of particles
flp.c = 0.01 * ones( geom.tnp, 1 );
% Compute the pressure as speed of sound times density
flp.p = (flp.c.^2).*flp.rho;      

end