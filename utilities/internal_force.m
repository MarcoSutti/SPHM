function [ flp, sph, int_dedt, int_dvdt ] = internal_force( geom, sph, flp, ia )

% function [ flp, sph, int_dedt, int_dvdt ] = internal_force( geom, sph, flp, ia )
% Calculates the internal forces on the RHS of the Navier-Stokes equations,
% i.e., the pressure gradient and the gradient of the viscous stress
% tensor, used by the time integration.
% Moreover, the entropy production due to viscous dissipation and the change
% of internal energy per mass, de/dt, are calculated.

% art_dvdt: part of the acceleration due to the internal forces

% Created:     ??.??.2011
% Last change: 10.06.2021

%   Jun 9, 2021:
%       Check and minor changes.
%   Aug 22, 2021:
%       V = geom.v';  % 22.08.2019: makes the vectorized for loop faster.


% Initialization of the deviator of the strecthing tensor, velocity divergence,
% viscous energy, internal energy, art_dvdt ...
% devD_xx(i), devD_yy(i), devD_zz(i), devD_xy(i), etc. are the components
% of the deviator of the strecthing tensor for particle i
devD_xx = zeros( geom.tnp, 1 );
devD_yy = zeros( geom.tnp, 1 );
devD_zz = zeros( geom.tnp, 1 );
devD_xy = zeros( geom.tnp, 1 );
devD_xz = zeros( geom.tnp, 1 );
devD_yz = zeros( geom.tnp, 1 );
sph.viscous_entropy = zeros( geom.tnp, 1 );
flp.c = zeros( geom.tnp, 1 );
p_work = zeros( geom.tnp, 1 );
% int_dedt = zeros( geom.tnp, 1 );       % time rate of change of the specific internal energy
int_dvdt = zeros( geom.tnp, geom.dim );  % acceleration due to internal forces

%==========================================================================
% Calculate the SPH approximation for the deviator of the rate of
% deformation tensor (aka stretching tensor):
%       devD_ab = (1/2) * (v_a,b + v_b,a) - (1/3) * vc,c * delta_ab.
if sph.visc                      % If the fluid is viscous, then...
    
    v_tr = geom.v';  % 22.08.2019: makes the vectorized for loop faster.
    
    for k=1:ia.niap   % ... for each interacting pair...
        i = ia.pair_i(k);          % ... take the index of the first particle in the pair...
        j = ia.pair_j(k);          % ... take the index of the second particle in the pair...
        
        % Calculate the velocity difference (aka relative velocity)
        v_ji = v_tr(:,j) - v_tr(:,i);
        
        % Compute the contribution of the current (k-th) interacting pair
        % to the components of devD.
        % NB: h_xx, h_xy, etc. are all dummy variables.
        if geom.dim == 1
            % one-dimensional devD has 1 component
            h_xx = 2/3 * v_ji(1) * ia.dW_ijdx(k,1);
            % ia.dW_ijdx is the matrix of size ia.niap * geom.dim which
            % stores the derivatives of the kernel in the d-direction for every interacting pair
        elseif geom.dim == 2
            % two-dimensional devD has 4 components
            % (but 2 of them are equal, so we are left with 3 components)
            h_xx = 2/3 * v_ji(1) * ia.dW_ijdx(k,1) -  1/3 * v_ji(2) * ia.dW_ijdx(k,2);
            h_xy = 1/2 * ( v_ji(1) * ia.dW_ijdx(k,2) + v_ji(2) * ia.dW_ijdx(k,1) );
            h_yy = 2/3 * v_ji(2) * ia.dW_ijdx(k,2) - 1/3 * v_ji(1) * ia.dW_ijdx(k,1);
        elseif geom.dim == 3
            % three-dimensional devD has 9 components
            % (but 6 of them are equal, so we are left with 6 components)
            h_xx = 2/3 * v_ji(1) * ia.dW_ijdx(k,1) ...
                - 1/3 * ( v_ji(2) * ia.dW_ijdx(k,2) + v_ji(3) * ia.dW_ijdx(k,3) );
            h_xy = 1/2 * ( v_ji(1) * ia.dW_ijdx(k,2) + v_ji(2) * ia.dW_ijdx(k,1) );
            h_xz = 1/2 * ( v_ji(1) * ia.dW_ijdx(k,3) + v_ji(3) * ia.dW_ijdx(k,1) );
            h_yy = 2/3 * v_ji(2) * ia.dW_ijdx(k,2) ...
                - 1/3 * ( v_ji(1) * ia.dW_ijdx(k,1) + v_ji(3) * ia.dW_ijdx(k,3) );
            h_yz = 1/2 * ( v_ji(2) * ia.dW_ijdx(k,3) + v_ji(3) * ia.dW_ijdx(k,2) );
            h_zz = 2/3 * v_ji(3) * ia.dW_ijdx(k,3) ...
                - 1/3 * ( v_ji(1) * ia.dW_ijdx(k,1) + v_ji(2) * ia.dW_ijdx(k,2) );
        end
        
        
        % Now for the current interacting pair, for both particles i and j
        % perform the summation to find all the components of devD
        
        % Precompute the multiplying factors:
        factor_i = sph.mass(i)/flp.rho(i);
        factor_j = sph.mass(j)/flp.rho(j);
        
        if geom.dim == 1
            devD_xx(i) = devD_xx(i) + factor_j * h_xx;
            devD_xx(j) = devD_xx(j) + factor_i * h_xx;
        elseif geom.dim == 2
            devD_xx(i) = devD_xx(i) + factor_j * h_xx;
            devD_xx(j) = devD_xx(j) + factor_i * h_xx;
            devD_xy(i) = devD_xy(i) + factor_j * h_xy;
            devD_xy(j) = devD_xy(j) + factor_i * h_xy;
            devD_yy(i) = devD_yy(i) + factor_j * h_yy;
            devD_yy(j) = devD_yy(j) + factor_i * h_yy;
        elseif geom.dim == 3
            devD_xx(i) = devD_xx(i) + factor_j * h_xx;
            devD_xx(j) = devD_xx(j) + factor_i * h_xx;
            devD_xy(i) = devD_xy(i) + factor_j * h_xy;
            devD_xy(j) = devD_xy(j) + factor_i * h_xy;
            devD_xz(i) = devD_xz(i) + factor_j * h_xz;
            devD_xz(j) = devD_xz(j) + factor_i * h_xz;
            devD_yy(i) = devD_yy(i) + factor_j * h_yy;
            devD_yy(j) = devD_yy(j) + factor_i * h_yy;
            devD_yz(i) = devD_yz(i) + factor_j * h_yz;
            devD_yz(j) = devD_yz(j) + factor_i * h_yz;
            devD_zz(i) = devD_zz(i) + factor_j * h_zz;
            devD_zz(j) = devD_zz(j) + factor_i * h_zz;
        end
        
        % MS, 10.09.2021: Is there something missing here?
%         hvcc = dot( dvx, dwdx );
%         vcc(i) = vcc(i) + mass(j) * hvcc/rho(j);
%         vcc(j) = vcc(j) + mass(i) * hvcc/rho(i);        
        
    end
end

%==========================================================================
% Calculate sph.viscous entropy, pressure and speed of sound for
% every particle
%==========================================================================
if sph.visc      % ... if we are dealing with a sph.viscous fluid, then...
    % ... we compute the sum of the squared components
    % of the deviator of D
    if geom.dim == 1
        sum_devD_comp_sq = devD_xx.^2;
    elseif geom.dim == 2
        sum_devD_comp_sq = devD_xx.^2 + 2 * devD_xy.^2 + devD_yy.^2;
        % Because of symmetry, devD_xy = devD_yx
    elseif geom.dim == 3
        sum_devD_comp_sq = devD_xx.^2 + devD_yy.^2 + devD_zz.^2 ...
            + 2 * ( devD_xy.^2 + devD_xz.^2 + devD_yz.^2 );
    end
    
    % Viscous entropy variation (p. 120):
    % sph.viscous_entropy = 2*sph.eta/flp.rho*(sum of devD components squared)
    % MS, 09.06.2021: OK
    sph.viscous_entropy = 2 * (sph.eta./flp.rho).*sum_devD_comp_sq;
end

% Compute pressure and speed of sound with the equation of state
if flp.fluid_type == 1           % If we are dealing with a gas
    flp = p_gas( flp, sph );
elseif flp.fluid_type == 2       % If we are dealing with water
    flp = p_art_water( flp, sph, geom );
end

%==========================================================================
% Calculate the SPH approximation for the equation of motion,
% i.e. for pressure force -p,a/rho
% and viscous force (eta * 2devD_ab),b/rho
% and the SPH approximation of the internal energy change de/dt due to the
% pressure work -p/rho vc,c
%==========================================================================
for k=1:ia.niap   % For each interacting pair...
    i = ia.pair_i(k);          % ... take the index of the first particle in the pair...
    j = ia.pair_j(k);          % ... take the index of the second particle in the pair...
    
    p_work_part = 0;        % initialization
    %======================================================================
    %... if we use the symmetrization with 1/( flp.rho(i) * flp.rho(j) )
    if sph.part_approx == 1
        
        rho_ij = 1/( flp.rho(i) * flp.rho(j) );
        
        for d=1:geom.dim     %... for every spatial geom.dimension...
            
            % eom_forces: forces appearing in the equation of motion
            eom_forces = -( flp.p(i) + flp.p(j) ) * ia.dW_ijdx(k,d);
            
            % p_work_part: divergence of the velocity difference times (p(i) + p(j))
            p_work_part = p_work_part + ( geom.v(j,d) - geom.v(i,d) ) * eom_forces;
            
            if sph.visc    % If the fluid is sph.viscous then...
                %... to the pressure forces we add the sph.viscous forces...
                if d==1          %... x-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xx(i) ...
                        + sph.eta(j) * devD_xx(j) ) * ia.dW_ijdx(k,1);
                    if geom.dim >= 2        % if the problem is 2D, we add the y-part to the x-coordinate of acceleration
                        eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xy(i) ...
                            + sph.eta(j) * devD_xy(j) ) * ia.dW_ijdx(k,2);
                        if geom.dim==3      % if the problem is 3D, we add the z-part to the x-coordinate of acceleration
                            eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xz(i) ...
                                + sph.eta(j) * devD_xz(j) ) * ia.dW_ijdx(k,3);
                        end
                    end
                elseif d==2      %... y-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xy(i) ...
                        + sph.eta(j) * devD_xy(j) ) * ia.dW_ijdx(k,1) ...
                        + 2 * ( sph.eta(i) * devD_yy(i) ...
                        + sph.eta(j) * devD_yy(j) ) * ia.dW_ijdx(k,2);
                    if geom.dim == 3    %if the problem is 3D, we add the z-part to the y-coordinate of acceleration
                        eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_yz(i) ...
                            + sph.eta(j) * devD_yz(j) ) * ia.dW_ijdx(k,3);
                    end
                elseif d==3      %... z-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xz(i) ...
                        + sph.eta(j) * devD_xz(j) ) * ia.dW_ijdx(k,1) ...
                        + 2 * ( sph.eta(i) * devD_yz(i) + sph.eta(j) * devD_yz(j) ) * ia.dW_ijdx(k,2) ...
                        + 2 * ( sph.eta(i) * devD_zz(i) + sph.eta(j) * devD_zz(j) ) * ia.dW_ijdx(k,3);
                end
            end
            
            %Acceleration due to internal forces (int_dvdt) in the d-direction...
            %... for particle i:
            int_dvdt(i,d) = int_dvdt(i,d) + sph.mass(j) * rho_ij * eom_forces;
            %... for particle j:
            int_dvdt(j,d) = int_dvdt(j,d) - sph.mass(i) * rho_ij * eom_forces;
            % minus sign because the forces on particle j are opposite to those on
            % particle i
            
        end       % end of the loop over d
        
        % Pressure work appearing in the energy equation...
        % ... for particle i:
        p_work(i) = p_work(i) + 0.5 * sph.mass(j) * rho_ij * p_work_part;
        % ... for particle j:
        p_work(j) = p_work(j) + 0.5 * sph.mass(i) * rho_ij * p_work_part;
        
        
        %======================================================================
        % ... if we use the symmetrization with 1/flp.rho(i)^2, 1/flp.rho(j)^2
    elseif sph.part_approx == 2
        for d=1:geom.dim
            eom_forces = -( flp.p(i)/flp.rho(i)^2 + flp.p(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,d);
            p_work_part = p_work_part + ( geom.v(j,d) - geom.v(i,d) ) * eom_forces;
            
            % Viscous force
            if sph.visc            %If the fluid is sph.viscous, then...
                %... to the pressure forces we add the sph.viscous forces...
                if d==1        %... x-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xx(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_xx(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,1);
                    if geom.dim >= 2      %if the problem is 2D, we add the y-part to the x-coordinate of acceleration
                        eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xy(i)/flp.rho(i)^2 ...
                            + sph.eta(j) * devD_xy(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,2);
                        if geom.dim == 3  %if the problem is 3D, we add the z-part to the x-coordinate of acceleration
                            eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xz(i)/flp.rho(i)^2 ...
                                + sph.eta(j) * devD_xz(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,3);
                        end
                    end
                    
                elseif d==2    %... y-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xy(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_xy(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,1) ...
                        + 2 * ( sph.eta(i) * devD_yy(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_yy(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,2);
                    if geom.dim == 3  %if the problem is 3D, we add the z-part to the y-coordinate of acceleration
                        eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_yz(i)/flp.rho(i)^2 ...
                            + sph.eta(j)*devD_yz(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,3);
                    end
                    
                elseif d==3    %... z-coordinate of acceleration...
                    eom_forces = eom_forces + 2 * ( sph.eta(i) * devD_xz(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_xz(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,1) ...
                        + 2 * ( sph.eta(i) * devD_yz(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_yz(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,2) ...
                        + 2 * ( sph.eta(i) * devD_zz(i)/flp.rho(i)^2 ...
                        + sph.eta(j) * devD_zz(j)/flp.rho(j)^2 ) * ia.dW_ijdx(k,3);
                end
            end
            
            % Momentum equation in the d-direction...
            % ... for particle i:
            int_dvdt(i,d) = int_dvdt(i,d) + sph.mass(j) * eom_forces;
            % ... for particle j:
            int_dvdt(j,d) = int_dvdt(j,d) - sph.mass(i) * eom_forces;
            % minus sign because the forces on particle j are opposite to
            % those on particle i
        end       % end of the loop over d
        
        % Pressure work appearing in the energy equation...
        %... for particle i:
        p_work(i) = p_work(i) + (1/2) * sph.mass(j) * p_work_part;
        %... for particle j:
        p_work(j) = p_work(j) + (1/2) * sph.mass(i) * p_work_part;
    end
    
end     %end of loop over the pairs


%==========================================================================
% Change of specific internal energy 
%==========================================================================
int_dedt = p_work + sph.viscous_entropy;

end