%==========================================================================
% Contains all the parameters for the SPH simulation.
% Created:     17.06.2021
% Last change: 20.06.2021

%   Jun 10, 2021:
%       Added option for quartic smoothing kernel by Liu, Liu and Lam
%       (2002).
%==========================================================================
%                              INPUT DATA
%==========================================================================
% geom.dim : Dimension of the problem (1, 2 or 3)
geom.dim = 1;
%==========================================================================
tip.nstart = 0;

% Time-step type (1 = fixed, 2 = adaptive, Monaghan [1994], 3 = adaptive,
% Monaghan and Kos [1999])
tip.ts_type = 1;

% If you choose an adaptive time-step technique, this value is used just
% for the first time-step and all the subsequent time-steps are computed
% according to the chosen adaptive technique.
% N.B.: the time-step has to satisfy the Courant-Friedrichs-Levy (CFL)
% condition for an explicit scheme to be stable.
tip.dt = 5e-3;

% Maximum number of time-steps
tip.max_nts = 100;
%==========================================================================
% Type of fluid (1 = gas, 2 = water)
flp.fluid_type = 1;
%==========================================================================
% plt.real_time = false, no real time
%               = true, real time
plt.real_time = true;

% Type of plot:
% plt.type = 1, position (and velocity quiver if plt.velocity = 1)
%          = 2, velocity colormap
%          = 3, pressure colormap
plt.type = 1;

% plt.velocity = false, no velocity quiver
%              = true, velocity quiver
plt.velocity = false;
%==========================================================================

% Maximum number of real particles
geom.max_nrp = 12000;

% SPH algorithm for particle approximation (sph.part_approx)
% sph.part_approx = 1 : (e.g. (p(i)+p(j))/(flp.rho(i)*flp.rho(j))
%                   2 : (e.g. (p(i)/flp.rho(i)^2+p(j)/flp.rho(j)^2)
sph.part_approx = 2;

% Nearest neighbor particle searching (nnps) method
% nnps = 1 : Simplest and direct searching
%        2 : Sorting grid linked list
%        3 : Tree algorithm
sph.NNPS = 1;

% Smoothing length evolution (sph.sle) algorithm
% sph.sle = 0 : Keep unchanged,
%           1 : h(i)=sigma*(sph.mass(i)/flp.rho(i))^(1/geom.dim)        [Monaghan, 2005]
%           2 : dh/tip.dt = (-1/geom.dim)*(h/flp.rho)*(drho/tip.dt)         [Benz, 1989]
sph.sle = 0;

% Factor for defining the smoothing length.
% Monaghan and Kajtar (2009) also used a value of 1.3-1.5
hsml_factor = 1.5;

% Smoothing kernel function
% sph.skf = 1, Cubic spline kernel by W4 - Spline (Monaghan 1985)
%         = 2, Gauss kernel   (Gingold and Monaghan 1981)
%         = 3, Quintic kernel (Morris 1997)
%         = 4, Quartic kernel (Liu, Liu and Lam 2002)
sph.skf = 4;

%--------------------------------------------------------------------------
% Switches for different senarios
%--------------------------------------------------------------------------
% summation_density = true : Use density summation model in the code,
%                     false: Use continuity equation
sph.sum_density = true;

% sph.average_velocity = true : Monaghan treatment on average velocity,
%                        false: No average treatment.
sph.avg_velocity = false;

% sph.virtual_part     = true : Use virtual particle,
%                        false: No use of vritual particle.
sph.virtual_part = false;

% sph.visc             = true : Consider sph.viscosity,
%                        false: No sph.viscosity.
sph.visc = false;

% sph.ext_force        = true : Consider external force,
%                        false: No external force.
sph.ext_force = false;    % gravity force and boundary particle force

% sph.visc_artificial   = true : Consider artificial sph.viscosity,
%                         false: No considering of artificial sph.viscosity.
sph.visc_artificial = true;

% sph.heat_artificial = true : Consider artificial heating,
%                       false: No considering of artificial heating.
sph.heat_artificial = false;

% sph.gravity_force = true : Considering sph.gravity_force,
%                     false: No considering of sph.gravity_force
sph.gravity_force = false;

% Gravitational acceleration
sph.g = -9.81;

% sph.normalized_density = true : Density normalization by using CSPM,
%                        = false: No normalization.
sph.normalized_density = false;

% sph.boundary_force_approach
% Opt. 1: Boundary particle force and penalty anti-penetration force
% [Monaghan, 1994] in Lennard-Jones form.
% Opt. 2: Normal Boundary Repulsive force according to Monaghan & Kos (1999).
% Opt. 3: Normal Boundary Repulsive force according to Monaghan (2003).
sph.boundary_force_approach = 1;

% MS, 10.06.2021:
% Maximum number of interactions:
sph.max_nia = 100 * geom.max_nrp;

% Verbose option:
% sph.verbose = 0, no information on console
%             = 1, basic information
%             = 2, Statistics: interactions per particle
sph.verbose = 1;

% Parameters for the artificial viscosity: alpha and beta are constants
% that are typically set around one (see p. 126).
% "The choice of alpha and beta is not critical, but they should be near
% alpha = 1 and beta = 2 for best results" [cited from Monaghan, 1992].
% In Liu and Liu's code they are both set equal to 1.
sph.alpha = 1;      % The viscosity associated with alpha produces a bulk viscosity
sph.beta = 1;       % Suppresses particle interpenetration at high Mach number
%==========================================================================
% Boundary particle force and penalty anti-penetration force [Monaghan, 1994]
% in Lennard-Jones form.
%==========================================================================
% Cutoff distance
% Pay attention in choosing the cutoff distance. If it is too large, some
% particles may feel the repulsive force from the boundary already in the
% initial distribution. If too small, real particles may penetrate the
% boundary before feeling the repulsive force.
% In light of the initial particle distribution we have assumed, the most
% appropriate value should be about one half of the particle spacing.
% sph.r0 = 0.48;

% MS, 15.06.2021: In Liu & Liu's code, this is set to 1.25e-5
sph.rr0 = 1.25e-5;
                
% D is a coefficient to be chosen considering the physical configuration.
% "D is a problem dependent parameter, and should be chosen to be in the
% same scale as the square of the largest velocity" [cited from Liu & Liu, 2003, p. 141]
% For problems involving dams, bores, weirs with fluid depth H, we take
% D = 5gH, but also D = 10gH or D = gH [cited from Monaghan, 1994]
% sph.D = abs( sph.g * ( geom.y_max - geom.y_min ) );
% sph.D = abs( sph.g * ( 1e-3 ) );
% sph.n_1 must be > sph.n_2
% (possible choices: sph.n_1 = 4, sph.n_2 = 2; sph.n_1 = 12, sph.n_2 = 6)
% In Liu and Liu's code: sph.n_1 = 12, sph.n_2 = 4.
sph.n_1 = 12;
sph.n_2 = 4;

% Factor appearing in the penalty anti-penetration force.
sph.dd = 1e-2;

%==========================================================================
%                           END OF INPUT DATA
%==========================================================================