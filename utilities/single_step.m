function [ geom, flp, tip ] = single_step( geom, sph, flp, tip )

% function [ geom, flp, tip ] = single_step( geom, sph, flp, tip )
% Determines the right hand side of a differential equation in a single
% step for performing time integration.
% In this routine and its subroutines the SPH algorithms are performed.
% Created:     ??.??.2011
% Last change: 24.06.2020

%   Jun 24, 2021:
%       Added "geom" as output parameter.
%       Added if statement for generation of virtual particles.
%   Jun 10, 2021:
%       Minor changes. Introduction of the array irp.

% Initializations.
art_visc = zeros( geom.nrp, 1 );
art_heat = zeros( geom.nrp, 1 );
ext_dvdt = zeros( geom.nrp, geom.dim );
art_dvdt = zeros( geom.nrp, geom.dim );

geom.nvp = 0;

% Position of virtual boundary particles.
if sph.virtual_part
    switch sph.example
        case 2
            [ geom, sph, flp ] = shear_cavity_virtual_part( geom, sph, flp );
        case 3
            [ geom, sph, flp ] = dam_collapse_virtual_part( geom, sph, flp );
        otherwise
            error('Unknown example.');
    end
end

% Find the nearest neighbouring particles, i.e. the particles falling
% into the support domain of a given particle. This process, called NNPS,
% is necessary at every time-step, since the particle distribution evolves
% with time.
if sph.NNPS == 1
    % Call the direct find/pairwise interaction algorithm
    ia = direct_find( geom, sph );
elseif sph.NNPS == 2
    % Call the linked-list algorithm
    error('Linked-list algorithm not yet implemented!')
elseif sph.NNPS == 3
    % Call the tree search algorithm
    error('Tree-search algorithm not yet implemented!')
end

% Particle approximation of density
if sph.sum_density
    % if we are using the summation density approach
    flp = summation_density( geom, sph, flp, ia );
else
    % if we are using the continuity density approach
    tip.drhodt = continuity_density( geom, sph, ia );
end

% Dynamic viscosity:
sph = viscosity( geom, sph );

% Internal forces:
[ flp, sph, int_dedt, int_dvdt ] = internal_force( geom, sph, flp, ia );

% Artificial viscosity:
if sph.visc_artificial
    [ art_dvdt, art_visc, max_phi_ij ] = artificial_viscosity( geom, sph, flp, ia );
end

% External forces:
if sph.ext_force
    [ ext_dvdt, magnitude_ext_dvdt ] = external_force( geom, sph, ia );
end

%==========================================================================
% Choose the technique for the time-step
if tip.ts_type == 2
    % Calculate the adaptive time-step according to Monaghan [1994].
    % Liu and Liu, p. 142.
    dt_cv = min( sph.hsml./( flp.c + 0.6 * ( sph.alpha * flp.c + sph.beta * max_phi_ij ) ) );
    dt_f = min( sqrt( sph.hsml./magnitude_ext_dvdt ) );
    tip.dt = min( 0.4 * dt_cv, 0.25 * dt_f );
elseif tip.ts_type == 3
    % Adaptive time step according to Monaghan and Kos [1999]:
    tip.dt = min( 0.3 * sph.hsml./( flp.c + max_phi_ij ) );
end
%==========================================================================

% Update the smoothing length
sph = hsml_update( geom, sph, flp, tip, ia );

% Calculating artificial heat
if sph.heat_artificial
    art_heat = artificial_heat( geom, sph, flp, ia );
end

% Calculating average velocity of each particle to avoid particle penetration
if sph.avg_velocity
    avg_v = avg_velocity( geom, sph, flp, ia );
    tip.avg_v = avg_v;
end


% NB: ONLY FOR REAL PARTICLES!
% Indexes of Real Particles:
irp = 1:geom.nrp;
% The TOTAL ACCELERATION is given by:
% - the acceleration due to internal forces, int_dvdt
% - the acceleration due to external forces, ext_dvdt
% - the acceleration due to artificial sph_param.viscosity, art_dvdt
tip.dvdt = int_dvdt(irp, :) + ext_dvdt(irp, :) + art_dvdt(irp, :);

% The TOTAL ENERGY VARIATION is given by:
% - energy variation due to internal forces, int_dedt
% - artificial sph_param.viscosity, art_visc
% - artificial heat, art_heat
tip.dedt = int_dedt(irp) + art_visc(irp) + art_heat;

end