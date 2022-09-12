%==========================================================================
% Generates the profile plots for the shock tube example.

% Created:     20.06.2021
% Last change: 15.07.2022

%   Jul 15, 2022:
%       Added save epsc.
%   Jun 20, 2021:
%       Created.
%==========================================================================
% Startup
SPHM_startup;
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
xlim = 0.4;

% Load simulation data from matfile:
fileName = 'results/example_1_max_nts_40.mat';
load( fileName );
%--------------------------------------------------------------------------
% Compute exact solution
st_param.x0 = -0.6;
st_param.x5 = 0.6;

st_param.rhoL = 1.0;
st_param.pL   = 1.0;
st_param.uL   = 0.0;

st_param.rhoR = 0.25;
st_param.pR   = 0.1795;
st_param.uR   = 0.0;

st_param.gamma = 1.4;

st_param.t = 0.20;

% Position of the membrane:
st_param.x_barrier = 0.0;

%--------------------------------------------------------------------------

st_param.N = 50;

st_results = shocktube_exact( st_param );


%--------------------------------------------------------------------------
% Plots
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');
% figure('units','normalized','outerposition',[0 0 1 1])

% Density profile.
figure(1);
plot( x_hist(:, end), rho_hist(:, end), 'Color', plt.color.green, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('$\rho$ [kg/m$^{3}$]')
axis( [ -xlim, xlim, 0, 1.2 ] )
grid on
hold on
plot( st_results.x_vec_1, st_results.rho1 * ones(length(st_results.x_vec_1),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_2, st_results.rho2, '--k', 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.rho3 * ones(length(st_results.x_vec_3),1), '--k', 'LineWidth', 2 )
plot( st_results.x3*[1, 1], [st_results.rho3, st_results.rho4], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.rho4 * ones(length(st_results.x_vec_4),1), '--k', 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.rho4, st_results.rho5], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.rho5 * ones(length(st_results.x_vec_5),1), '--k', 'LineWidth', 2 )

% Legend
legend( {'SPH solution', 'exact solution'}, ...
    'FontSize', 11, 'Location', 'northeast' );


% Save plot to eps file
fileName = 'plots/density_profile';
saveas( gcf, fileName, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName);


% Pressure profile.
figure(2);
plot( x_hist(:, end), p_hist(:, end), 'Color', plt.color.blue, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('$p$ [N/m$^{2}$]')
axis( [ -xlim, xlim, 0, 1.2 ] )
grid on
hold on
plot( st_results.x_vec_1, st_results.p1 * ones(length(st_results.x_vec_1),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_2, st_results.p2, '--k', 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.p3 * ones(length(st_results.x_vec_3),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.p4 * ones(length(st_results.x_vec_4),1), '--k', 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.p4, st_results.p5], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.p5 * ones(length(st_results.x_vec_5),1), '--k', 'LineWidth', 2 )

% Legend
legend( {'SPH solution', 'exact solution'}, ...
    'FontSize', 11, 'Location', 'northeast' );

% Save plot to eps file
fileName = 'plots/pressure_profile';
saveas( gcf, fileName, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName);


% Velocity profile.
figure(3);
plot( x_hist(:, end), v_hist(:, end), 'Color', plt.color.red, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('$v$ [m/s]')
axis( [ -xlim, xlim, -0.025, 1.05 ] )
grid on
hold on
plot( st_results.x_vec_1, st_results.u1 * ones(length(st_results.x_vec_1),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_2, st_results.u2, '--k', 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.u3 * ones(length(st_results.x_vec_3),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.u4 * ones(length(st_results.x_vec_4),1), '--k', 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.u4, st_results.u5], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.u5 * ones(length(st_results.x_vec_5),1), '--k', 'LineWidth', 2 )

% Legend
legend( {'SPH solution', 'exact solution'}, ...
    'FontSize', 11, 'Location', 'northeast' );

% Save plot to eps file
fileName = 'plots/velocity_profile';
saveas( gcf, fileName, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName);


% Internal energy profile.
figure(4);
plot( x_hist(:, end), e_hist(:, end), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('$e$ [J/kg]')
axis( [ -xlim, xlim, 1.6, 2.75 ] )
grid on
hold on
plot( st_results.x_vec_1, st_results.e1 * ones(length(st_results.x_vec_1),1), '--k', 'LineWidth', 2 )
plot( st_results.x_vec_2, st_results.e2, '--k', 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.e3 * ones(length(st_results.x_vec_3),1), '--k', 'LineWidth', 2 )
plot( st_results.x3*[1, 1], [st_results.e3, st_results.e4], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.e4 * ones(length(st_results.x_vec_4),1), '--k', 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.e4, st_results.e5], '--k', 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.e5 * ones(length(st_results.x_vec_5),1), '--k', 'LineWidth', 2 )

% Legend
legend( {'SPH solution', 'exact solution'}, ...
    'FontSize', 11, 'Location', 'northeast' );

% Save plot to eps file
fileName = 'plots/energy_profile';
saveas( gcf, fileName, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName);
