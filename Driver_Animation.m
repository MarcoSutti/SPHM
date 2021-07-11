%==========================================================================
% This script can be used to generate the animation corresponding to the
% time evolution of SPH particles.

% Created:     20.06.2021
% Last change: 09.07.2021

%   Jul 9, 2021:
%       Added commands to save the animation in a video file.
%   Jun 20, 2021:
%       Created.
%==========================================================================
% Startup
sph_startup;
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Type of plot:
% plt.type = 1, position (and velocity quiver if plt.velocity = 1)
%          = 2, velocity colormap
%          = 3, pressure colormap
plt.type = 1;

% plt.velocity = 0, no velocity quiver
%              = 1, velocity quiver
plt.velocity = 1;

stride = 1;

example = 1;
max_nts = 40;

%--------------------------------------------------------------------------

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Load data                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName = ['results/example_', num2str(example), '_max_nts_', ...
    num2str(max_nts), '.mat'];

% Load simulation data from matfile:
load( fileName );


fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Animation                          |\n');
fprintf('+--------------------------------------------------------------+\n');

figure('units','normalized','outerposition',[0 0 1 1])

videoFileName = ['videos/example_', num2str(example), '_max_nts_', ...
    num2str(max_nts), '_stride_', num2str(stride), '.avi'];

v = VideoWriter( videoFileName, 'Motion JPEG AVI' );
v.Quality = 95;
open(v)

% Do animation:
for its = 1:stride:size(x_hist,2)/geom.dim
    
    if geom.dim == 1
        if plt.type == 1
            % Plot the particle position
            plot( x_hist(:, its), zeros(geom.nrp, 1), 'o', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceColor', ...
                plt.color.Sky, 'MarkerSize', 7 );
            
            if plt.velocity
                hold on
                % Plot the velocity field
                quiver( x_hist(:, its), zeros(geom.nrp, 1), ...
                    v_hist(:, its), zeros(geom.nrp, 1), ...
                    1.5, 'LineWidth', 1.25, 'Color',  plt.color.blue )
            end
        elseif plt.type == 2
            scatter( x_hist(:, its), zeros(geom.nrp, 1), 50, ...
                abs(v_hist(:, its)), 'filled', 'MarkerEdgeColor',  plt.color.Sky )
            % colorbar;
        elseif plt.type == 3
            scatter( x_hist(:, its), zeros(geom.nrp, 1), 50, ...
                p_hist(:, its), 'filled', 'MarkerEdgeColor',  plt.color.Sky )
            % colorbar;
        end
        axis equal
        axis( [ -2.5, 1, -0.5, 0.5 ] )
        title( [ 'time = ', num2str(t_hist(its)) ] );
        
        pause(0.01)
        hold off
    else
        if plt.type == 1
            % Plot the particle position
            handle_rp = plot( x_hist(:, 2*(its-1)+1 ), x_hist(:, 2*its ), ...
                'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',  plt.color.Sky, 'MarkerSize', 7 );
            
            if plt.velocity
                hold on
                % Plot the velocity field
                quiver( x_hist(:, 2*(its-1)+1 ), x_hist(:, 2*its ), ...
                    v_hist(:, 2*(its-1)+1 ), v_hist(:, 2*its ), 1.5, ...
                    'LineWidth', 1.25, 'Color',  plt.color.blue )
            end
        elseif plt.type == 2
            
            speed = sqrt( v_hist(:, 2*(its-1)+1 ).^2 + v_hist(:, 2*its ).^2 );
            
            scatter( x_hist(:, 2*(its-1)+1 ), x_hist(:, 2*its ), ...
                50, speed, 'filled', 'MarkerEdgeColor',  plt.color.Sky )
            % colorbar;
        elseif plt.type == 3
            scatter( x_hist(:, 2*(its-1)+1 ), x_hist(:, 2*its ), 50, ...
                p_hist(:, its), 'filled', 'MarkerEdgeColor',  plt.color.Sky )
            % colorbar;
        end
        axis equal
        axis( [ -0.05e-3, 1.05e-3, -0.05e-3, 1.05e-3 ] )
        %         axis( [ geom.x_min - geom.wf_left, geom.x_max + geom.wf_right, ...
        %             geom.y_min - geom.wf_down, geom.y_max + geom.wf_up ] )
        hold on
        
        plot( geom.x_bp, geom.y_bp, 's', 'MarkerEdgeColor',  plt.color.gray3, ...
            'MarkerFaceColor',  plt.color.gray1, 'MarkerSize', 7 );
        
        title( [ 'time = ', num2str(its) ] );
        
        pause(0.01)
        hold off
    end
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v)