%==========================================================================
% Generates the profile plots for the shock tube example.

% Created:     20.06.2021
% Last change: 20.06.2021

%   Jun 20, 2021:
%       Created.
%==========================================================================
% Startup
sph_startup;
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');

% Load simulation data from matfile:
fileName = 'results/example_1_max_nts_40.mat';
load( fileName );

figure('units','normalized','outerposition',[0 0 1 1])

% Density profile.
subplot(2,2,1);
plot( x_hist(:, end), rho_hist(:, end), 'Color', plt.color.green, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('density $\rho$ [kg/m$^{3}$]')
axis( [ -0.4, 0.4, 0, 1.2 ] )
grid on


% Pressure profile.
subplot(2,2,2);
plot( x_hist(:, end), p_hist(:, end), 'Color', plt.color.blue, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('pressure $p$ [N/m$^{2}$]')
axis( [ -0.4, 0.4, 0, 1.2 ] )
grid on


% Velocity profile.
subplot(2,2,3);
plot( x_hist(:, end), v_hist(:, end), 'Color', plt.color.red, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('velocity $v$ [m/s]')
axis( [ -0.4, 0.4, -0.05, 1 ] )
grid on


% Internal energy profile.
subplot(2,2,4);
plot( x_hist(:, end), e_hist(:, end), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('internal energy $e$ [J/kg]')
% axis( [ -0.4, 0.4, 1.8, 2.6 ] )
grid on
