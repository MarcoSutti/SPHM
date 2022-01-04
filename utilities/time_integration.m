function [ geom, sph, flp  ] = time_integration( geom, sph, flp, tip, plt )

% function [ geom, sph, flp  ] = time_integration( geom, sph, flp, tip, plt )
% Purpose: Performs time integration.

% Created:     ??.??.2011
% Last change: 24.06.2021

%   Jun 24, 2021:
%       Added "geom" as output parameter of single_step.
%   Jun 20, 2021:
%       Added save_output.
%   Jun 15, 2021:
%       Added temp_e and temp_rho.
%   Jun 10, 2021:
%       Introduction of the array irp.

if sph.verbose > 0
    fprintf('--------------------------------------------------------------\n');
    fprintf(' Starting SPH time integration\n');
    fprintf('--------------------------------------------------------------\n');
end

% Initializations
time = tip.nstart;
current_ts = tip.nstart;
dedt = zeros( geom.nrp, 1 );
e_min = zeros( geom.nrp, 1 );
rho_min = zeros( geom.nrp, 1 );
v_min = zeros( geom.nrp, geom.dim );

% Indexes of Real Particles:
irp = 1:geom.nrp;

% Initialization of arrays where we store the history of position,
% velocity, pressure, energy, density, and time.
x_hist = zeros( geom.nrp, geom.dim * tip.max_nts );
v_hist = zeros( geom.nrp, geom.dim * tip.max_nts );
p_hist = zeros( geom.nrp, tip.max_nts );
e_hist = zeros( geom.nrp, tip.max_nts );
rho_hist = zeros( geom.nrp, tip.max_nts );
t_hist = zeros( tip.max_nts, 1 );

for its = ( tip.nstart + 1 ):( tip.nstart + tip.max_nts )
    
    current_ts = current_ts + 1;
    
    if sph.verbose > 0
        fprintf(' timestep: %5d\n', current_ts );
    end
    
    if its ~= 1
        %If this is not the first time step, then set minimum energy
        e_min = sph.e(irp);
        
        temp_e = zeros( geom.nrp, 1 );
        
        if geom.dim == 1
            temp_e = -geom.nsym * flp.p.*geom.v(irp)./( geom.x(irp).*flp.rho );
        end
        
        %... update energy half a time step
        sph.e(irp) = sph.e(irp) + (tip.dt/2) * ( dedt + temp_e );
        
        %If the energy so computed is negative, then set it to zero
        sph.e( sph.e(irp) < 0 ) = 0;
        
        if ~sph.sum_density
            %If we are not using the summation density approach...
            %... set minimum density
            rho_min = flp.rho(irp);
            
            temp_rho = zeros( geom.nrp, 1 );
            
            if geom.dim == 1
                temp_rho = -geom.nsym * flp.rho.*geom.v(irp, 1)./geom.x(irp, 1);
            end
            
            %... update density half a time step
            flp.rho(irp) = flp.rho(irp) + (tip.dt/2) * ( tip.drhodt(irp) + temp_rho );
        end
        
        %... set minimum velocity...
        v_min = geom.v(irp, :);
        %... update velocity half a time step
        geom.v(irp, :) = geom.v(irp, :) + (tip.dt/2) * tip.dvdt(irp, :);
    end
    
    % Definition of variables out of the function vector:
    % Call the function single_step
    [ geom, flp, tip ] = single_step( geom, sph, flp, tip );
    
    if its == 1
        temp_e = zeros( geom.nrp, 1 );
        
        if geom.dim == 1
            temp_e = -geom.nsym * flp.p.*geom.v(irp)./( geom.x(irp).*flp.rho );
        end
        
        % If this is the first timestep, then update energy half a time step
        sph.e(irp) = sph.e(irp) + (tip.dt/2) * ( dedt + temp_e );
        
        % If the energy so computed is negative, then set it to zero
        sph.e( sph.e(irp) < 0 ) = 0;
        
        if ~sph.sum_density
            % If we are not using the summation density approach...
            
            temp_rho = zeros( geom.nrp, 1 );
            
            if geom.dim == 1
                temp_rho = -geom.nsym * flp.rho.*geom.v(irp, 1)./geom.x(irp, 1);
            end
            
            % ... update density half a time step
            flp.rho(irp) = flp.rho(irp) + (tip.dt/2) * ( tip.drhodt(irp) + temp_rho );
        end
        
        % Update velocity half a time step
        if sph.avg_velocity == true
            geom.v(irp, :) = geom.v(irp, :) + (tip.dt/2) * tip.dvdt(irp, :) ...
                + tip.avg_v(irp, :);
        elseif sph.avg_velocity == false
            geom.v(irp, :) = geom.v(irp, :) + (tip.dt/2) * tip.dvdt(irp, :);
        end
        
        % Update particle position a full time step
        if geom.dim ~= 1
            geom.x(irp, :) = geom.x(irp, :) + tip.dt * geom.v(irp, :);
        else
            geom.x(irp) = geom.x(irp) + tip.dt * geom.v(irp);
        end
    else
        % Else, if this is not the first time step, then update energy a
        % full time step
        temp_e = zeros( geom.nrp, 1 );
        
        if geom.dim == 1
            temp_e = -geom.nsym * flp.p.*geom.v(irp)./( geom.x(irp).*flp.rho );
        end
        
        sph.e(irp) = e_min + tip.dt * ( dedt + temp_e );
        
        % If the energy so computed is negative then set it to zero
        sph.e( sph.e(irp) < 0 ) = 0;
        
        if ~sph.sum_density
            % If we are not using the summation density approach...
            
            temp_rho = zeros( geom.nrp, 1 );
            
            if geom.dim == 1
                temp_rho = -geom.nsym * flp.rho.*geom.v(irp, 1)./geom.x(irp, 1);
            end
            
            % If we are not using the summation density approach update
            % density a full time step
            flp.rho(irp) = rho_min + tip.dt * ( tip.drhodt(irp) + temp_rho );
        end
        
        % update velocity a full time step...
        if sph.avg_velocity == true
            geom.v(irp, :) = v_min + tip.dt * tip.dvdt(irp, :) + tip.avg_v(irp, :);
        elseif sph.avg_velocity == false
            geom.v(irp, :) = v_min + tip.dt * tip.dvdt(irp, :);
        end
        
        % update particle position a full time step...
        geom.x(irp, :) = geom.x(irp, :) + tip.dt * geom.v(irp, :);
    end
    
    % calculate current time
    time = time + tip.dt;
    
    % MS, 20.06.2021: Save position, velocity, pressure, energy and density
    % histories to the variables x_hist, v_hist, p_hist, e_hist, and rho_hist.
    if geom.dim ~= 1
        x_hist( :, 2*(its-1)+1:2*its ) = geom.x(irp, :);
        v_hist( :, 2*(its-1)+1:2*its ) = geom.v(irp, :);
        p_hist( :, its ) = flp.p(irp, :);
        e_hist( :, its ) = sph.e(irp, :);
        rho_hist( :, its ) = flp.rho(irp, :);
    else
        x_hist( :, its ) = geom.x(irp);
        v_hist( :, its ) = geom.v(irp);
        p_hist( :, its ) = flp.p(irp);
        e_hist( :, its ) = sph.e(irp);
        rho_hist( :, its ) = flp.rho(irp);
    end
    
    t_hist( its, 1 ) = time;
    
    if plt.real_time
        plot_particle_evolution( geom, plt, time );
    end
    
end

if sph.verbose > 0
    fprintf('--------------------------------------------------------------\n');
    fprintf(' End of SPH time integration\n');
    fprintf('--------------------------------------------------------------\n');
end

tip.nstart = current_ts;

% MS, 20.06.2021: Save simulation data for later.
save_output( x_hist, v_hist, p_hist, e_hist, rho_hist, t_hist, geom, sph, tip );

end