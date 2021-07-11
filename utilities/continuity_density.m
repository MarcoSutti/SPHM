function [ drhodt ] = continuity_density( geom, sph, ia )

% function [ drhodt ] = continuity_density( geom, sph, ia )
% Calculates the density with SPH continuity approach.
% Created:     10.01.2021
% Last change: 10.06.2021

%   Jun 10, 2021:
%       Minor changes. Vectorization of the computation of W_i and
%       fluid_param.rho.
%       Use MATLAB's arrayfun to apply the function kernel to each element
%       of array sph.hsml.

% Inizialization of drhodt
drhodt = zeros( geom.tnp, 1 );

v_tr = geom.v';  % 22.08.2019: makes the vectorized for loop faster.

for k=1:ia.niap
    i = ia.pair_i(k);
    j = ia.pair_j(k);
    
    % Calculate the velocity difference (aka relative velocity)
    v_ji = v_tr(:,j) - v_tr(:,i);
    
    vcc = dot( v_ji, ia.dW_ijdx(k,:) );
    
    drhodt(i) = drhodt(i) + sph.mass(j) * vcc;
    drhodt(j) = drhodt(j) + sph.mass(i) * vcc;
end

end