function [ geom, sph, flp ] = dam_collapse_virtual_part( geom, sph, flp )

% function [ geom, sph, flp ] = dam_collapse_virtual_part( geom, sph, flp )
% Purpose: Generates the virtual boundary particles for the dam collapse
%          problem.

% Created:     06.08.2021
% Last change: 06.08.2021

% Virtual particle spacing
geom.vps = geom.rps/4;

%Initialization of the number of virtual particles
geom.nvp = 75/geom.vps+27.5/geom.vps+2.5/geom.vps;

% Store total number of particles
geom.tnp = geom.nrp + geom.nvp;

%==========================================================================
% Create geometry of domain:
%==========================================================================
% LOWER BOUNDARY
%==========================================================================
stride1 = 1:75/geom.vps;
geom.x(stride1 + geom.nrp, 1) = geom.vps * (stride1-1);
geom.x(stride1 + geom.nrp, 2) = 0;
geom.v(stride1 + geom.nrp, 1) = 0;
geom.v(stride1 + geom.nrp, 2) = 0;

%==========================================================================
% LEFT BOUNDARY
%==========================================================================
stride2 = 0:27.5/geom.vps;
geom.x(75/geom.vps + stride2 + geom.nrp, 1) = 0;
geom.x(75/geom.vps + stride2 + geom.nrp, 2) = geom.vps * stride2;
geom.v(75/geom.vps + stride2 + geom.nrp, 1) = 0;
geom.v(75/geom.vps + stride2 + geom.nrp, 2) = 0;

%==========================================================================
% TILTED BAR
%==========================================================================
stride3 = 0:2.5/geom.vps;
geom.x(75/geom.vps + 27.5/geom.vps + stride3 + geom.nrp, 1) = 50 + geom.vps * stride3;
geom.x(75/geom.vps + 27.5/geom.vps + stride3 + geom.nrp, 2) = geom.vps * stride3;
geom.v(75/geom.vps + 27.5/geom.vps + stride3 + geom.nrp, 1) = 0;
geom.v(75/geom.vps + 27.5/geom.vps + stride3 + geom.nrp, 2) = 0;

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