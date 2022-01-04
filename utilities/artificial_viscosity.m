function [ art_dvdt, art_visc, max_phi_ij ] = artificial_viscosity( geom, sph, flp, ia )

% function [ art_dvdt, art_visc, max_phi_ij ] = artificial_viscosity( geom, sph, flp, ia )
% art_dvdt: part of the acceleration due to the artificial viscosity

% Created:     ??.??.2011
% Last change: 09.06.2021

%   Jun 9, 2021:
%       Added some comments.
%   Aug 23, 2019:
%       There is a more efficient way to do the vectorization, see comments
%       in the code.

% Initialization of art_dvdt and art_visc
art_dvdt = zeros( geom.tnp, geom.dim );
art_visc = zeros( geom.tnp, 1 );

%==========================================================================
% Calculate SPH sum for artificial sph.viscosity
%==========================================================================
phi = zeros( ia.niap, 1 );

% MS, 22.08.2019: Why is the vectorized version slower than the original
% for loop? Answer: There is a wrong and a correct way to do the
% vectorization!!!
% It is more efficient by doing the operations on the TRANSPOSED array:
if geom.dim ~= 1
    v_tr = geom.v';
    x_tr = geom.x';
end

for k = 1:ia.niap          % For every interacting pair...
    i = ia.pair_i(k);
    j = ia.pair_j(k);
    
    if geom.dim == 1
        x_ij = geom.x(i) - geom.x(j);
        v_ij = geom.v(i) - geom.v(j);
    else
        % Position difference
        x_ij = x_tr(:,i) - x_tr(:,j);
        % Velocity difference (aka relative velocity)
        v_ij = v_tr(:,i) - v_tr(:,j);
    end
    
    % Calculate scalar product between v_ij and x_ij
    v_ij_dot_x_ij = dot( v_ij, x_ij );

    hsml_ij = ( sph.hsml(i) + sph.hsml(j) )/2;
    
    phi_ij = ( hsml_ij * v_ij_dot_x_ij )/( norm(x_ij)^2 + 0.01 * hsml_ij^2 );
    
    % Store all the phi_ij for a later use (in the calculation of the
    % adaptive time-step)
    phi(k) = phi_ij;
    
    %======================================================================
    % Artificial sph.viscous force only if v_ij * r_ij < 0
    %======================================================================
    if v_ij_dot_x_ij < 0
        % If the scalar product between v_ij and x_ij is negative (it means
        % compression), then calculate
        %       Pi_ij = ( -alpha*c_ij*phi_ij + beta*phi_ij^2 )/rho_ij
        % according to [Lattanzio et al., 1986; Monaghan, 1989]
                    
        c_ij   = ( flp.c(i) + flp.c(j) )/2;       % arithmetic mean of the speeds of sound
        rho_ij = ( flp.rho(i) + flp.rho(j) )/2;   % arithmetic mean of the densities
    
        % Monaghan type artificial viscosity (p. 126)
        Pi_ij  = ( - sph.alpha * c_ij * phi_ij + sph.beta * phi_ij^2)/rho_ij;
        
        %==================================================================
        % Add SPH approximation for artificial sph.viscous force to
        % the pressure terms in SPH equations
        % Equation of Motion:
        %... subtract the term to the acceleration of particle i
        art_dvdt(i,:) = art_dvdt(i,:) - Pi_ij * sph.mass(j) * ia.dW_ijdx(k,:);
        
        % ... and add it to the acceleration of particle j in order
        % to satisfy Newton's Third Law
        art_dvdt(j,:) = art_dvdt(j,:) + Pi_ij * sph.mass(i) * ia.dW_ijdx(k,:);
        
        %Artificial sph.viscosity (to be added to the energy equation):
        v_ij_dot_dW_ijdx = dot( v_ij, ia.dW_ijdx(k,:) );
        
        art_visc(i) = art_visc(i) + 0.5 * Pi_ij * sph.mass(j) * v_ij_dot_dW_ijdx;
        art_visc(j) = art_visc(j) + 0.5 * Pi_ij * sph.mass(i) * v_ij_dot_dW_ijdx;
    end
end

% Stores the maximum phi_ij (to be used in the calculation of the adaptive
% time-stepping)
max_phi_ij = max(phi);

end