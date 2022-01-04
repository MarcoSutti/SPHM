function [ sph ] = viscosity( geom, sph )

% function [ sph ] = viscosity( geom, sph )
% Purpose: Defines the fluid particle sph.viscosity eta [kg/(m*s)]

% Created:     ??.??.2011
% Last change: 09.06.2021

%   Jun 9, 2021:
%       The logical indexing should give the same result as the older loop.

% Initialization
sph.eta = zeros( geom.tnp, 1 );

% For gas particles, we do not have viscosity.
% For water particles, we assume the value of 1.0e-3 kg/(m*s)
sph.eta( abs(geom.part_type) == 2 ) = 1.0e-3;

end