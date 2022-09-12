function [ W, dWdx ] = kernel( r_ij, x_ij, hsml_ij, geom, sph )

% function [ W, dWdx ] = kernel( r_ij, x_ij, hsml_ij, geom, sph )
% Calculates the smoothing kernel W and its derivatives dWdx.

% Created:     2011
% Last change: 22.06.2021

%   Jun 22, 2021:
%       Added Liu, Liu and Lam quartic smoothing function.
%   Jun 9, 2021:
%       Removed W and dWdx set to zero at the end of each if statement,
%       because they are already initialized as zero arrays.

R = r_ij/hsml_ij;   % relative distance between particles i and j (see p. 62)

% Initialization of W and its derivatives:
W = 0;
dWdx = zeros( 1, geom.dim );

%==========================================================================
%      Cubic spline kernel by Monaghan and Lattanzio (1985) (see p. 64)
%==========================================================================
if sph.skf == 1
    if geom.dim == 1
        factor = 1/hsml_ij;
    elseif geom.dim == 2
        factor = 15/(7*pi*hsml_ij^2);
    elseif geom.dim == 3
        factor = 3/(2*pi*hsml_ij^3);
    end
    if R>=0 && R<=1
        W = factor * (2/3 - R^2 + 0.5 * R^3);
        % Evaluation of the derivative
        dWdx = factor * (-2 + 3/2*R)/hsml_ij^2 * x_ij;
    elseif R>1 && R<=2
        W = factor * 1/6 * (2-R)^3;
        % Evaluation of the derivative
        dWdx = -factor * 1/2 * ((2-R)^2)/hsml_ij * (x_ij/r_ij);
    end
%==========================================================================
%        Gauss kernel by Gingold and Monaghan (1981) (see p. 62-63)
%==========================================================================
elseif sph.skf == 2
    factor = 1/(pi^(geom.dim/2)*hsml_ij^geom.dim);
    if R>=0 && R<=3
        W = factor * exp(-R^2);
        dWdx = W * (-2*x_ij/hsml_ij^2);
    end
%==========================================================================
%               Quintic kernel by Morris (1997) (see p. 65-66)
%==========================================================================
elseif sph.skf == 3
    if geom.dim == 1
        factor = 1/(120*hsml_ij);
    elseif geom.dim == 2
        factor = 7/(478*pi*hsml_ij^2);
    elseif geom.dim == 3
        factor = 3/(359*pi*hsml_ij^3);
    end 
    if R>=0 && R<=1
        W = factor * ( (3-R)^5 - 6*(2-R)^5 + 15*(1-R)^5 );
        dWdx = factor * (-120 + 120*R - 50*R^2 )/hsml_ij^2 * x_ij;
    elseif R>=1 && R<=2
        W = factor * ( (3-R)^5 -6*(2-R)^5 );
        dWdx = factor * ( -5*(3-R)^4 + 30*(2-R)^4 )/hsml_ij * (x_ij/r_ij);
    elseif R>2 && R<=3
        W = factor * (3-R)^5;
        dWdx = factor * ( -5*(3-R)^4 )/hsml_ij * (x_ij/r_ij);
    end
%==========================================================================
%           Quartic kernel by Liu, Liu and Lam (2002) (see p. 92)
%==========================================================================
elseif sph.skf == 4
    if geom.dim == 1
        factor = 1/hsml_ij;
    elseif geom.dim == 2
        factor = 15/(7*pi*hsml_ij^2);
    elseif geom.dim == 3
        factor = 315/(208*pi*hsml_ij^3);
    end
    if R>=0 && R<=2
        W = factor * (2/3 - (9/8) * R^2 + (19/24) * R^3 - (5/32) * R^4);
        dWdx = factor * (-9/4 + (19/8) * R - (5/8) * R^2)/hsml_ij^2 * x_ij;
    end
else
    error('Please enter a valid choice of smoothing kernel.')
end
end