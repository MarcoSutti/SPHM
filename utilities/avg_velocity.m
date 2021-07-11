function [ avg_v ] = avg_velocity( geom, sph, flp, ia )

% Calculates the average velocity to correct velocity in order to prevent
% penetration (Monaghan, 1992). See p. 138, about XSPH technique,
% or Monaghan, 1992.

% avg_v: average velocity of each particle.
  
% epsilon --- a small constant chosen by experience, may lead to instability.
% 0<=epsilon<=1.0
% In most circumstances, epsilon=0.3 seems to be a good choice.
% In the simulation of explosions also epsilon>=0.5 can be used.

epsilon = 0.3;

%Initialization of the average velocity
avg_v = zeros( geom.tnp, geom.dim );   % MS, 24.06.2021: Shouldn't it be geom.nrp?

%Compute the average velocity
for k = 1:ia.niap         %For each interacting pair...
    i = ia.pair_i(k);
    j = ia.pair_j(k);
    
    % Arithmetic mean of the density
    rho_ij = ( flp.rho(i) + flp.rho(j) )/2;
    
    % Velocity difference (aka relative velocity)
    v_ij = geom.v(i,:) - geom.v(j,:);
    
    % Average velocity for particle i
    avg_v(i,:) = avg_v(i,:) - sph.mass(j)/rho_ij * v_ij * ia.W_ij(k);
    
    % Average velocity for particle j
    avg_v(j,:) = avg_v(j,:) + sph.mass(i)/rho_ij * v_ij * ia.W_ij(k);
end

avg_v = epsilon * avg_v;

end