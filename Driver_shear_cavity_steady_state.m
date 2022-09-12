%==========================================================================
% This script can be used to generate the plot of the shear-driven cavity
% problem at the steady state with particle positions and velocity quivers.

% Created:     2022.09.12
% Last change: 2022.09.12

%   Sep 12, 2022:
%       Created.
%==========================================================================
% Startup
SPHM_startup;
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Type of plot:
% plt.type = 1, position
%          = 2, velocity colormap
%          = 3, pressure colormap
plt.type = 1;

% plt.velocity = 0, no velocity quiver
%              = 1, velocity quiver
plt.velocity = 1;

example = 2;
max_nts = 10000;

time_idx = max_nts;

xlim = 1e-3;

%--------------------------------------------------------------------------

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Load data                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName = ['results/example_', num2str(example), '_max_nts_', ...
    num2str(max_nts), '.mat'];

% Load simulation data from matfile:
load( fileName );


fprintf('+--------------------------------------------------------------+\n');
fprintf('|                              Plot                            |\n');
fprintf('+--------------------------------------------------------------+\n');

figure('units','normalized','outerposition',[0 0 1 1]);

if plt.type == 1
    if plt.velocity == 0
        % Plot the particle position
        handle_rp = plot( x_hist(:, time_idx-1 ), x_hist(:, time_idx ), ...
            'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor',  plt.color.Sky, 'MarkerSize', 7 );
        
    elseif plt.velocity == 1
        % hold on
        % Plot the velocity field
        quiver( x_hist(:, time_idx-1 ), x_hist(:, time_idx ), ...
            v_hist(:, time_idx-1 ), v_hist(:, time_idx ), 1.5, ...
            'LineWidth', 1.25, 'Color',  plt.color.blue )
    end
elseif plt.type == 2
    
    speed = sqrt( v_hist(:, time_idx-1 ).^2 + v_hist(:, time_idx ).^2 );
    
    scatter( x_hist(:, time_idx-1 ), x_hist(:, time_idx ), ...
        50, speed, 'filled', 'MarkerEdgeColor',  plt.color.Sky )
    % colorbar;
elseif plt.type == 3
    scatter( x_hist(:, time_idx-1 ), x_hist(:, time_idx ), 50, ...
        p_hist(:, its), 'filled', 'MarkerEdgeColor',  plt.color.Sky )
    % colorbar;
end
axis equal
hold on

plot( geom.x_bp, geom.y_bp, 's', 'MarkerEdgeColor',  plt.color.gray3, ...
    'MarkerFaceColor',  plt.color.gray1, 'MarkerSize', 7 );

title_string = [ 'time = ', num2str(t_hist(time_idx)), ' s' ];
title( title_string );

axis( [ 0, xlim, 0, xlim ] )

xlabel('$ x_{1} $')
ylabel('$ x_{2} $')

% axis off

% Save plot to eps file
fileName = 'plots/shear_cavity_steady_state';
saveas( gcf, fileName, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName);
