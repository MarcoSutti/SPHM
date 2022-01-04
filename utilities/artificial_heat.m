function [ art_heat ] = artificial_heat( geom, sph, flp, ia )

% function [ art_heat ] = artificial_heat( geom, sph, flp, ia )
% Computes the artificial heat [Noh, 1987; Fulk, 1994; Monaghan, 1995].
% The artificial heat conduction term computed here is added to the energy
% equation.

% Parameter for the artificial heat conduction:
alpha = 0.1;
beta = 1.0;

% Initializations.
div_v = zeros( geom.tnp, 1 );
art_heat = zeros( geom.tnp, 1 );

%==========================================================================
%SPH approximation of the velocity divergence
for k = 1:ia.niap   % For every interacting pair...
    i = ia.pair_i(k);
    j = ia.pair_j(k);
    
    % Velocity difference for the current pair
    v_ji = geom.v(j,:) - geom.v(i,:);
    
    % Divergence of the velocity difference for the current pair
    div_v_ji = dot( v_ji, ia.dW_ijdx(k,:) );
    
    % SPH approximation of the velocity divergence at particle i
    div_v(i) = div_v(i) + sph.mass(j)/flp.rho(j) * div_v_ji;
    
    % SPH approximation of the velocity divergence at particle j
    div_v(j) = div_v(j) + sph.mass(i)/flp.rho(i) * div_v_ji;
end

for k = 1:ia.niap        % For each interacting pair...
    i = ia.pair_i(k);    % ... take the index of the first particle in the pair...
    j = ia.pair_j(k);    % ... take the index of the second particle in the pair...
    
    hsml_ij= ( sph.hsml(i) + sph.hsml(j) )/2;    % arithmetic mean of the smoothing lengths
    rho_ij = ( flp.rho(i) + flp.rho(j) )/2;      % arithmetic mean of the densities
    
    % Position difference:
    x_ij = geom.x(i,:) - geom.x(j,:);
    
    % Norm of x_ij squared:
    r_ij_sq = norm( x_ij )^2;
    
    % Scalar product between x_ij and ia.dW_ijdx(k,:)
    x_ij_dot_dW_ijdx  = dot( x_ij, ia.dW_ijdx(k,:) );
    
    q_i = alpha * sph.hsml(i) * flp.rho(i) * flp.c(i) * abs( div_v(i) ) ...
        + beta * sph.hsml(i)^2 * flp.rho(i) * div_v(i)^2;
    q_j = alpha * sph.hsml(j) * flp.rho(j) * flp.c(j) * abs( div_v(j) ) ...
        + beta * sph.hsml(j)^2 * flp.rho(j) * div_v(j)^2;
    q_ij = ( q_i + q_j )/2;
    h = 2 * q_ij/rho_ij * 1/( r_ij_sq + 0.01 * hsml_ij^2 ) * x_ij_dot_dW_ijdx;   % dummy variable
    
    %Artificial heat (to be added to the energy equation)
    art_heat(i) = art_heat(i) + h * ( sph.e(i) - sph.e(j) );
    art_heat(j) = art_heat(j) + h * ( sph.e(j) - sph.e(i) );
    
end
end