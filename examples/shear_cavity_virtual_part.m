function [ geom, sph, flp ] = shear_cavity_virtual_part( geom, sph, flp )

% function [ geom, sph, flp ] = shear_cavity_virtual_part( geom, sph, flp )
% Purpose: Generates the virtual boundary particles for the 2D shear cavity
% driven problem.

% Created:     22.06.2021
% Last change: 06.08.2021

%   Jan 4, 2022:
%       Added the variables ivp_###.

% nrp_per_side = 40;
xl = 1e-3;
v_inf = 1e-3;

nvp_per_side = 2*sqrt(geom.nrp);

%Initialization of the number of virtual particles
geom.nvp = 4 * nvp_per_side;

% Store total number of particles
geom.tnp = geom.nrp + geom.nvp;

% Virtual particle spacing
geom.vps = geom.rps/2;

%==========================================================================
% Create geometry of domain:
%==========================================================================
% UPPER BOUNDARY
%==========================================================================
stride1 = 1:nvp_per_side;
ivp_upper = stride1 + geom.nrp;
geom.x(ivp_upper, 1) = geom.vps * (stride1-1);
geom.x(ivp_upper, 2) = xl;
geom.v(ivp_upper, 1) = v_inf;
geom.v(ivp_upper, 2) = 0;

%==========================================================================
% LOWER BOUNDARY
%==========================================================================
stride2 = nvp_per_side+1:2*nvp_per_side;
ivp_lower = stride2 + geom.nrp;
geom.x(ivp_lower, 1) = geom.vps * (stride1-1);
geom.x(ivp_lower, 2) = 0;
geom.v(ivp_lower, 1) = 0;
geom.v(ivp_lower, 2) = 0;

%==========================================================================
% LEFT BOUNDARY
%==========================================================================
stride3 = 2*nvp_per_side+1:3*nvp_per_side;
ivp_left = stride3 + geom.nrp;
geom.x(ivp_left, 1) = 0;
geom.x(ivp_left, 2) = geom.vps * stride1;
geom.v(ivp_left, 1) = 0;
geom.v(ivp_left, 2) = 0;

%==========================================================================
% RIGHT BOUNDARY
%==========================================================================
stride4 = 3*nvp_per_side+1:geom.nvp;
ivp_right = stride4 + geom.nrp;
geom.x(ivp_right, 1) = xl;
geom.x(ivp_right, 2) = geom.vps * stride1;
geom.v(ivp_right, 1) = 0;
geom.v(ivp_right, 2) = 0;

%==========================================================================
%         Store the physical variables for the virtual particles
%==========================================================================
% Indexes for Virtual Particles:
ivp = geom.nrp+1:geom.tnp;
flp.rho(ivp) = 1000;         % water density [kg/m^3]
sph.mass(ivp) = flp.rho(ivp) * geom.rps^2;
flp.p(ivp) = 0;
sph.e(ivp) = 357.1;

% Particle label to distinguish between fluid particles and boundary
% particles. We assign a positive label to real fluid particles and a
% negative label to boundary particles.
geom.part_type(ivp) = -flp.fluid_type;
sph.hsml(ivp) = geom.rps;

% Store virtual particles in x_bp and y_bp arrays:
geom.x_bp = geom.x(ivp, 1);
geom.y_bp = geom.x(ivp, 2);

end