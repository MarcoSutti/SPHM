function [ sph ] = hsml_update( geom, sph, flp, timeint, ia )

% function [ sph ] = hsml_update( geom, sph, flp, timeint, ia )
% Purpose: Calculates the evolution of the smoothing length.

% Created:     ??.??.2011
% Last change: 09.06.2021

%   Jun 9, 2021:
%       Check and minor changes.

if sph.sle == 0
    % Keep smoothing length unchanged.
    return
elseif sph.sle == 1
    % Evolve smoothing length according to [Monaghan, 2005]:
    % h(i) = sigma*(sph.mass(i)/flp.rho(i))^(1/geom.dim)
    sigma = 1.3;   % MS, 10.07.2021: In Liu & Liu's code, this is set to 2.
    sph.hsml = sigma * ( sph.mass./flp.rho ).^( 1/geom.dim );
elseif sph.sle == 2
    % Evolve smoothing length according to [Benz, 1989]:
    % dh/timeint.dt = (-1/geom.dim)*(h/rho)*(drho/timeint.dt)
    % Let's plug into it the continuity equation: drho/timeint.dt = -rho*div(v)
    % and we get: dh/timeint.dt = (h/geom.dim)*div(v)
    % where div(v) can be discretized using the SPH approximation.
    
    % Initialization of the velocity divergence
    div_v = zeros( geom.nrp, 1 );
    
    %================================================================
    % We construct the SPH approximation of the velocity divergence
    for k = 1:ia.niap         % For every interacting pair...
        i = ia.pair_i(k);     % ... take the index of the first particle in the pair...
        j = ia.pair_j(k);     % ... take the index of the second particle in the pair...
        
        % Velocity difference for the current pair
        v_ji = geom.v(j,:) - geom.v(i,:);
        
        % Divergence of the velocity difference for the current pair
        div_v_ji = dot( v_ji, ia.dW_ijdx(k,:) );
        
        % SPH approximation of the velocity divergence at particle i
        div_v(i) = div_v(i) + sph.mass(j)/flp.rho(j) * div_v_ji;
        
        % SPH approximation of the velocity divergence at particle j
        div_v(j) = div_v(j) + sph.mass(i)/flp.rho(i) * div_v_ji;
    end
    
    %================================================================
    % Calculate the time rate of change of the smoothing length
    dhsml_dt = (sph.hsml./geom.dim).*div_v;
    
    % Update the smoothing length
    sph.hsml = sph.hsml + timeint.dt * dhsml_dt;
    
    % If the updated smoothing length is negative, subtract timeint.dt*dhsml(i) instead of adding it
    % Use logical indexing
    idx_neg = sph.hsml<=0;
    sph.hsml(idx_neg) = sph.hsml(idx_neg) - timeint.dt * dhsml_dt(idx_neg);
    
end
end