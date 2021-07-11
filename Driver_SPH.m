%==========================================================================
% Driver SPH.
% Main script for running a SPH simulation.
% This code is based on the equivalent FORTRAN code of Liu GR, Liu MB.
% The main reference for this code is the book:
% [1] Liu GR, Liu MB (2003) Smoothed particle hydrodynamics: a meshfree
% particle method. World Scientific, Singapore.

% Created:     ??.??.2011
% Last change: 08.07.2021

%   Jul 8, 2021:
%       Added more comments.
%   Jun 10, 2021:
%       Added maximum number of interactions sph.max_nia and sph.verbose.
%==========================================================================
% Startup
sph_startup;
%==========================================================================
%                              INPUT DATA
%==========================================================================
% Choose an example.
sph.example = 2;
%==========================================================================
%                      GENERATION OF PROBLEM GEOMETRY
%==========================================================================
switch sph.example
    case 1
        shock_tube_parameters;
        [ geom, sph, flp ] = shock_tube( sph, plt );
    case 2
        shear_cavity_parameters;
        [ geom, sph, flp ] = shear_cavity( geom, sph, flp, plt );
    otherwise
        error('Unknown example.');
end
%==========================================================================
%                            TIME INTEGRATION
%==========================================================================
% Perform the time integration
[ geom, sph, flp ] = time_integration( geom, sph, flp, tip, plt );
