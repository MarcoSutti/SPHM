function [ ia ] = direct_find( geom, sph )

% function [ ia ] = direct_find( geom, sph )
% Implements the Nearest Neighboring Particle Searching (NNPS).
% This function makes use of a pairwise interaction technique (Hockney
% and Eastwood, 1988; Hernquist and Katz, 1989; Riffert et al., 1995) in
% order to reduce the computational effort (see p. 148).

% Created:     ??.??.2011
% Last change: 22.06.2021

%   Jun 22, 2021:
%       Replaced the if conditional statement with switch, case, otherwise.
%   Jun 10, 2021:
%       Added maximum number of interactions.

%==========================================================================
% Parameter for choosing the k constant appearing in the radius of the
% support of the smoothing function
% if sph.skf == 1
%     scale_k = 2;
% elseif sph.skf == 2
%     scale_k = 3;
% elseif sph.skf == 3
%     scale_k = 3;
% elseif sph.skf == 4
%     scale_k = 2;
% end

switch sph.skf
    case 1
        scale_k = 2;
    case 2
        scale_k = 3;
    case 3
        scale_k = 3;
    case 4
        scale_k = 2;
    otherwise
        error('Please select a valid value for the smoothing kernel function skf.')
end
%==========================================================================

% Initialization of nnp: number of neighboring particles for the
% i-th particle
nnp = zeros( geom.tnp, 1 );

if geom.dim ~= 1
    x_tr = geom.x';   % makes the vectorized loop faster to execute.
end

ia.niap = 0;       % initialization of n_interact_pairs: total number of interacting pairs
for i=1:(geom.tnp-1)
    for j=i+1:geom.tnp
        % Difference between x_i and x_j
        
        if geom.dim == 1
            x_ij = geom.x(i) - geom.x(j);
        else
            x_ij = x_tr(:,i) - x_tr(:,j);
        end
        
        % Distance between particle i and particle j
        r_ij = norm(x_ij);
        
        % Arithmetic mean of the smoothing lengths (simmetrization needed
        % if h varies in both space and time)
        hsml_ij = (1/2) * ( sph.hsml(i) + sph.hsml(j) );
        
        if r_ij < scale_k * hsml_ij                        % If particle j falls into the support of particle i,
            if ia.niap < sph.max_nia
                ia.niap = ia.niap + 1;      % then ij is counted as an interacting pair.
                ia.pair_i( ia.niap ) = i;               % stores the index of the first particle in the current interacting pair
                ia.pair_j( ia.niap ) = j;               % stores the index of the second particle in the current interacting pair
                nnp(i) = nnp(i) + 1;        % updates the number of neighboring particles for particle i
                nnp(j) = nnp(j) + 1;        % updates the number of neighboring particles for particle j
                
                % Calculate the smoothing function and its derivative in the
                % d-direction for the current interacting pair
                [ W, dWdx ] = kernel( r_ij, x_ij, hsml_ij, geom, sph );
                % Remember that dWdw is a vector of length d.
                
                % We store the derivative of the kernel function in
                % the d-direction for the current interacting pair
                ia.W_ij(ia.niap) = W;
                % interactions.W_ij is a vector of length ia.niap,
                % storing the value of the kernel for each interacting pair
                
                ia.dW_ijdx(ia.niap, :) = dWdx(:);
                % interactions.dW_ijdx is a matrix of size ia.niap times geom.dim,
                % which stores the derivatives of the kernel in the direction d
                % for each interacting pair
                
            else
                error('Too many interactions.')
            end
        end
    end
end

sum_iap = 0;
max_iap = 0;
min_iap = 1000;
no_iap = 0;

for i=1:geom.tnp
    sum_iap = sum_iap + nnp(i);
    if nnp(i) >= max_iap
        max_iap = nnp(i);
        max_p = i;
    end
    if nnp(i) <= min_iap
        min_iap = nnp(i);
        min_p = i;
    end
    if nnp(i) == 0
        no_iap = no_iap + 1;
    end
end

if sph.verbose > 1
    fprintf('Statistics: interactions per particle:\n')
    fprintf(' Particle: \t  %5d, %5d\n', max_p, max_iap );
    fprintf(' Particle: \t  %5d, %5d\n', min_p, min_iap );
    fprintf(' Total number of pairs: \t  %5d\n', ia.niap );
    fprintf(' Particles with no interactions: \t %5d\n', no_iap );
end

end