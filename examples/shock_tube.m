function [ geom, sph, flp ] = shock_tube( sph, plt )

% function [ geom, sph, flp ] = shock_tube( geom, sph, flp, plt )
% Generates the initial data for the 1D shock tube problem.
% See Example 3.6, p. 94, in Liu and Liu, 2003.

% Created:     17.06.2021
% Last change: 20.06.2021

%   Jun 17, 2021:
%       Created.

tube_length = 0.6;
geom.dim = 1;

% Symmetry of the problem
% geom.nsym = 0 : no symmetry,
%           = 1 : axis symmetry,
%           = 2 : center symmetry.
geom.nsym = 0;
geom.nrp = 400;
geom.nvp = 0;
geom.tnp = geom.nrp + geom.nvp;
geom.part_spac = tube_length/80;
flp.fluid_type = 1;   % ideal gas

sph.mass = 0.75/geom.nrp * ones( geom.nrp, 1 );
sph.hsml = 2 * geom.part_spac * ones( geom.nrp, 1 );
geom.part_type = flp.fluid_type * ones( geom.nrp, 1 );

% Define inital position and velocity:
geom.x = zeros( geom.nrp, 1 );
geom.v = zeros( geom.nrp, 1 );
sph.e = zeros( geom.nrp, 1 );
flp.rho = zeros( geom.nrp, 1 );
flp.p = zeros( geom.nrp, 1 );

% Particles in the high-density region:
geom.x(1:320) = -tube_length + geom.part_spac/4.*(0:319);
sph.e( geom.x <= 1e-8 ) = 2.5;
flp.rho( geom.x <= 1e-8 ) = 1;
flp.p( geom.x <= 1e-8 ) = 1;

% Particles in the low-density region:
geom.x(321:geom.nrp) = geom.part_spac.*(1:80);
sph.e( geom.x > 1e-8 ) = 1.795;
flp.rho( geom.x > 1e-8 ) = 0.25;
flp.p( geom.x > 1e-8 ) = 0.1795;

if plt.real_time
    % Plot.
    figure('units','normalized','outerposition',[0 0 1 1])
    plot( geom.x, 0*geom.x, 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', plt.color.Sky, 'MarkerSize', 7 );
    pause(0.05)
end

end