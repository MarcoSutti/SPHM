function [ ext_dvdt, magnitude_ext_dvdt ] = external_force( geom, sph, ia )

% function [ ext_dvdt, magnitude_ext_dvdt ] = external_force( geom, sph, ia )
% Calculates the external forces:
%     (1) gravitational forces;
%     (2) boundary forces exerted by virtual particles type I.

% ext_dvdt: part of the acceleration due to the external forces.

% Created:     ??.??.2011
% Last change: 24.06.2021

%   Jul 10, 2021:
%       Return magnitude_ext_dvdt.
%   Jun 24, 2021:
%       Restored original formulation at lines  34-61 according to Liu and
%       Liu's code.
%   Jun 15, 2021:
%       Added factor 0.01 in line 52.
%   Jun 9, 2021:
%       Check and minor changes.

% Initialization of ext_dvdt
ext_dvdt = zeros( geom.tnp, geom.dim );
magnitude_ext_dvdt = zeros( geom.tnp, 1 );

%==========================================================================
% Assign the gravitational force
%==========================================================================
if sph.gravity_force 
    % If we consider gravity force, then for each particle we assign g to
    % the acceleration in the geom.dim-direction
    ext_dvdt( :, geom.dim ) = sph.g;
end

%==========================================================================
% Assign the boundary forces
%==========================================================================
if sph.boundary_force_approach == 1
    %======================================================================
    % Opt. 1: Boundary particle force and penalty anti-penetration force
    % [Monaghan, 1994] in Lennard-Jones form.
    %======================================================================
    % r_0 = sph.r0 * geom.part_spac;
    
    for k=1:ia.niap             % For every interacting pair...
        i = ia.pair_i(k);       % ... take the index of the first particle in the pair...
        j = ia.pair_j(k);       % ... take the index of the second particle in the pair...
        
        % maybe this loop can also be vectorized...
        % i's = ia.pair_i(1:ia.niap)
        % j's = ia.pair_j(1:ia.niap)        
        
        % If particles i and j are of "different kind"...
        if geom.part_type(i) * geom.part_type(j) < 0
            % Compute the distance between i and j
            x_ij = geom.x(i,:) - geom.x(j,:);
            r_ij = norm(x_ij);
            
            % If the distance between i and j is less than the cutoff distance...
            if r_ij < sph.rr0
                %... compute the repulsive force per unit mass
                % f = sph.D * ( (r_0/r_ij)^sph.n_1 - (r_0/r_ij)^sph.n_2 )/r_ij^2;
                f = ( (sph.rr0/r_ij)^sph.n_1 - (sph.rr0/r_ij)^sph.n_2 )/r_ij^2;
                
                % Add the repulsive force to the acceleration
                ext_dvdt(i,:) = ext_dvdt(i,:) + sph.dd * f * x_ij;
            end
        end
    end
elseif sph.boundary_force_approach == 2
    %======================================================================
    % Opt. 2: Normal Boundary Repulsive force according to Monaghan & Kos (1999).
    %======================================================================
    for k=1:ia.niap              %For every interacting pair...
        i = ia.pair_i(k);                    %... take the index of the first particle in the pair...
        j = ia.pair_j(k);                    %... take the index of the second particle in the pair...
        if geom.part_type(i) * geom.part_type(j)<0     %... if i and j are of "different kind"...
            %disp('Interaction with the boundary detected')
            x_ij = geom.x(i,:) - geom.x(j,:); %x_ij(d): distance between i and j in the d-direction
            r_ij = norm(x_ij);            %distance between i and j
            
            if geom.part_type(i)<0    % If particle i is a boundary particle...
                index_bp = i;
                index_fp = j;
            else
                index_bp = j;
                index_fp = i;
            end
            
            % Mapping the index for all the particles to the index just for the
            % boundary particles:
            index_bp = index_bp - geom.tnp;
            
            % Calculation of the local coordinates x and y
            y_local = abs( dot( x_ij, geom.unit_normal(index_bp,:) ) );
            x_local = sqrt( r_ij^2 - y_local^2 );
            
            % Calculation of the function B(x,y)
            A = 1/sph.hsml(index_bp)*0.01*3500;   %c(index_fp)^2
            q = y_local/( 2*sph.hsml(index_bp) );
            
            % Calculate R(y)
            if q<1
                R = A/sqrt(q)*(1-q);
            else
                R = 0;
            end
            
            % Calculate P(x)
            if x_local < 3*geom.virt_part_spac
                P = 0.5*( 1 + cos(pi*x_local/geom.virt_part_spac) );
            else
                P = 0;
            end
            
            % Compute B(x,y) = R(y)*P(x):
            B = R*P;
            
            f = sph.mass(1)/( sph.mass(1) + sph.mass(end) )*B/r_ij;
            
            %Add the repulsive force to the acceleration in the d-direction
            ext_dvdt(i,:) = ext_dvdt(i,:) + f * x_ij;
        end
    end
elseif sph.boundary_force_approach == 3
    %======================================================================
    % Opt. 3: Normal Boundary Repulsive force according to Monaghan (2003).
    %======================================================================
    for k=1:ia.niap              %For every interacting pair...
        i = ia.pair_i(k);                    %... take the index of the first particle in the pair...
        j = ia.pair_j(k);                    %... take the index of the second particle in the pair...
        if geom.part_type(i) * geom.part_type(j) < 0    %... if i and j are of "different kind"...
            %disp('Interaction with the boundary detected')
            x_ij = geom.x(i,:) - geom.x(j,:); %x_ij(d): distance between i and j in the d-direction
            r_ij = norm(x_ij);            %distance between i and j
            
            if geom.part_type(i)<0    % If particle i is a boundary particle...
                index_bp = i;
                index_fp = j;
            else
                index_bp = j;
                index_fp = i;
            end
            
            % Mapping the index for all the particles to the index just for the
            % boundary particles:
            index_bp = index_bp - geom.nrp;
            
            % Calculation of the local coordinates x and y
            y_local = abs( dot( x_ij, geom.unit_normal(index_bp,:) ) );
            x_local = sqrt( r_ij^2 - y_local^2 );
            
            % Calculation of the function B(x,y)
            beta = 0.02 * 3500/y_local;    %c(index_fp)^2
            q = y_local/sph.hsml(index_bp);
            
            % Calculate Gamma(y)
            if q >= 0 && q <= 0.66
                Gamma = 0.66 * beta;
            elseif q > 0.66 && q < 1
                Gamma = beta * (2*q - 1.5*q^2);
            elseif q >= 1 && q <= 2
                Gamma = 0.5 * beta * (2-q)^2;
            else
                Gamma = 0;
            end
            
            % Calculate Chi(x)
            if x_local < geom.virt_part_spac
                Chi = 1 - x_local/geom.virt_part_spac;
            else
                Chi = 0;
            end
            
            % Compute B(x,y) = Gamma(y)*Chi(x):
            B = Gamma*Chi;
            
            f = sph.mass(1)/( sph.mass(1)+sph.mass(end) )*B/r_ij;
            
            % Add the repulsive force to the acceleration
            ext_dvdt(i,:) = ext_dvdt(i,:) + f * x_ij;
        end
    end
end

%==========================================================================
% Compute the modulus of ext_dvdt (i.e., the magnitude of the acceleration
% vector )
for i=1:geom.tnp
    magnitude_ext_dvdt(i) = norm( ext_dvdt(i,:) );
end

% Take the smallest magnitude of acceleration vector among all the particles
% (this value will be used in the computation of the adaptive time-step)
% min_modulus_ext_dvdt = min( magnitude_ext_dvdt );

end