function [ flp ] = summation_density( geom, sph, flp, ia )

% function [ flp ] = summation_density( geom, sph, flp, ia )
% Calculates the density with SPH summation algorithm.
% Created:     ??.??.2011
% Last change: 10.06.2021

%   Jun 10, 2021:
%       Vectorization of the computation of W_i and flp.rho.
%       Use MATLAB's arrayfun to apply the function kernel to each element
%       of array sph.hsml.

% Inizialization of hv   (dummy variable)
hv = zeros( geom.dim, 1 );

% Self density of each particle: W_i(i) (kernel for distance r=0)
r = 0;
% sd = zeros( geom.tnp, 1 );
% for i=1:geom.tnp
%     [ selfdens, ~ ] = kernel( r, hv, sph.hsml(i), geom, sph );
%     sd(i) = selfdens;
% end

% MS, 10.06.2021: Use MATLAB's arrayfun to apply the function kernel to
% each element of array sph.hsml
sd = arrayfun( @(hsml) kernel( r, hv, hsml, geom, sph ), sph.hsml );

% contribution of particle i to its own W_i
% (this can be looked on as an initialization of the kernel W_i(i))
W_i = sd.*sph.mass./flp.rho;

% Now to the self-contribution we add all the other contributions due to the
% interacting particles, to get the SPH approximation for W_i
for k=1:ia.niap   %for the k-th interacting pair...
    i = ia.pair_i(k);         %take the index of the first particle in the k-th interacting pair
    j = ia.pair_j(k);         %take the index of the first particle in the k-th interacting pair
    %The kernel of the first particle in the k-th interacting pair is
    %given by the self-contribution to the kernel plus the contribution
    %of the other particle in the pair, weighted with the kernel of the pair
    W_i(i) = W_i(i) + sph.mass(j)/flp.rho(j) * ia.W_ij(k);
    %Mutatis mutandis...
    W_i(j) = W_i(j) + sph.mass(i)/flp.rho(i) * ia.W_ij(k);
end

% Self-contribution to the density flp.rho
flp.rho = sd.*sph.mass;

% SPH summation density flp.rho
for k=1:ia.niap
    i = ia.pair_i(k);
    j = ia.pair_j(k);
    flp.rho(i) = flp.rho(i) + sph.mass(j) * ia.W_ij(k);
    flp.rho(j) = flp.rho(j) + sph.mass(i) * ia.W_ij(k);
end

% Calculate the normalized flp.rho [Randles and Libersky, 1996]
if sph.normalized_density
    flp.rho = flp.rho./W_i;
end

end