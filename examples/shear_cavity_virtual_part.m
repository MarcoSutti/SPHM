function [ geom, sph, flp ] = shear_cavity_virtual_part( geom, sph, flp, plt )

% function [ geom, sph, flp ] = shear_cavity_virtual_part( geom, sph, flp, plt )
% Purpose: Generates the virtual boundary particles for the 2D shear cavity
% driven problem.

% Created:     22.06.2021
% Last change: 24.06.2021

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
geom.x(stride1 + geom.nrp, 1) = geom.vps * (stride1-1);
geom.x(stride1 + geom.nrp, 2) = xl;
geom.v(stride1 + geom.nrp, 1) = v_inf;
geom.v(stride1 + geom.nrp, 2) = 0;

%==========================================================================
% LOWER BOUNDARY
%==========================================================================
stride2 = nvp_per_side+1:2*nvp_per_side;
geom.x(stride2 + geom.nrp, 1) = geom.vps * (stride1-1);
geom.x(stride2 + geom.nrp, 2) = 0;
geom.v(stride2 + geom.nrp, 1) = 0;
geom.v(stride2 + geom.nrp, 2) = 0;

%==========================================================================
% LEFT BOUNDARY
%==========================================================================
stride3 = 2*nvp_per_side+1:3*nvp_per_side;
geom.x(stride3 + geom.nrp, 1) = 0;
geom.x(stride3 + geom.nrp, 2) = geom.vps * stride1;
geom.v(stride3 + geom.nrp, 1) = 0;
geom.v(stride3 + geom.nrp, 2) = 0;

%==========================================================================
% RIGHT BOUNDARY
%==========================================================================
stride4 = 3*nvp_per_side+1:geom.nvp;
geom.x(stride4 + geom.nrp, 1) = xl;
geom.x(stride4 + geom.nrp, 2) = geom.vps * stride1;
geom.v(stride4 + geom.nrp, 1) = 0;
geom.v(stride4 + geom.nrp, 2) = 0;

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

% %==========================================================================
% %                       Plot the virtual particles
% %==========================================================================
% if plt.real_time
%     % Plot virtual particles
%     % figure( 'units', 'normalized', 'outerposition', [0 0 1 1] )
%     plot( geom.x_bp, geom.y_bp, 's', 'MarkerEdgeColor', plt.color.gray3, ...
%         'MarkerFaceColor', plt.color.gray1, 'MarkerSize', 7 );
%     axis equal
%     hold on
% 
% %     % Plot the unit normals to the boundaries
% %     quiver( geom.x_bp, geom.y_bp, geom.unit_normal(:,1), geom.unit_normal(:,2), ...
% %         1.5, 'LineWidth', 1.25, 'Color', plt.color.yellow )
%         
%     if plt.velocity
%         hold on
%         % Plot the velocity field
%         quiver( geom.x(:,1), geom.x(:,2), geom.v(:,1), geom.v(:,2), ...
%             1.5, 'LineWidth', 1.25, 'Color', plt.color.blue )
%     end
%     
%     title('Generation of boundary virtual particles')
%     pause
% end

end