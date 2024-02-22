function plot_shocktube_exact( st_results, st_param, plt )
%--------------------------------------------------------------------------
% PLOTS
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');

% figure('units','normalized','outerposition',[0 0 1 1])
figure(1)
% Density profile.
% subplot(2,2,1);
plot( st_results.x_vec_1, st_results.rho1 * ones(length(st_results.x_vec_1),1), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
hold on
plot( st_results.x_vec_2, st_results.rho2, 'Color', plt.color.Cerulean, 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.rho3 * ones(length(st_results.x_vec_3),1), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
plot( st_results.x3*[1, 1], [st_results.rho3, st_results.rho4], 'Color', plt.color.Cerulean, 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.rho4 * ones(length(st_results.x_vec_4),1), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.rho4, st_results.rho5], 'Color', plt.color.Cerulean, 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.rho5 * ones(length(st_results.x_vec_5),1), 'Color', plt.color.Cerulean, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('density $\rho$ [kg/m$^{3}$]')
axis( [ st_param.x0, st_param.x5, 0, 1.1 ] )
grid on


figure(2)
% Pressure profile.
plot( st_results.x_vec_1, st_results.p1 * ones(length(st_results.x_vec_1),1), 'Color', plt.color.red, 'LineWidth', 2 )
hold on
plot( st_results.x_vec_2, st_results.p2, 'Color', plt.color.red, 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.p3 * ones(length(st_results.x_vec_3),1), 'Color', plt.color.red, 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.p4 * ones(length(st_results.x_vec_4),1), 'Color', plt.color.red, 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.p4, st_results.p5], 'Color', plt.color.red, 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.p5 * ones(length(st_results.x_vec_5),1), 'Color', plt.color.red, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('pressure $p$ [N/m$^{2}$]')
axis( [ st_param.x0, st_param.x5, 0, 1.1 ] )
grid on


figure(3)
% Velocity profile.
plot( st_results.x_vec_1, st_results.u1 * ones(length(st_results.x_vec_1),1), 'Color', plt.color.orange, 'LineWidth', 2 )
hold on
plot( st_results.x_vec_2, st_results.u2, 'Color', plt.color.orange, 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.u3 * ones(length(st_results.x_vec_3),1), 'Color', plt.color.orange, 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.u4 * ones(length(st_results.x_vec_4),1), 'Color', plt.color.orange, 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.u4, st_results.u5], 'Color', plt.color.orange, 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.u5 * ones(length(st_results.x_vec_5),1), 'Color', plt.color.orange, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('velocity $v$ [m/s]')
axis( [ st_param.x0, st_param.x5, 0, 0.8 ] )
grid on


figure(4)
% Energy profile.
plot( st_results.x_vec_1, st_results.e1 * ones(length(st_results.x_vec_1),1), 'Color', plt.color.green, 'LineWidth', 2 )
hold on
plot( st_results.x_vec_2, st_results.e2, 'Color', plt.color.green, 'LineWidth', 2 )
plot( st_results.x_vec_3, st_results.e3 * ones(length(st_results.x_vec_3),1), 'Color', plt.color.green, 'LineWidth', 2 )
plot( st_results.x3*[1, 1], [st_results.e3, st_results.e4], 'Color', plt.color.green, 'LineWidth', 2 )
plot( st_results.x_vec_4, st_results.e4 * ones(length(st_results.x_vec_4),1), 'Color', plt.color.green, 'LineWidth', 2 )
plot( st_results.x4*[1, 1], [st_results.e4, st_results.e5], 'Color', plt.color.green, 'LineWidth', 2 )
plot( st_results.x_vec_5, st_results.e5 * ones(length(st_results.x_vec_5),1), 'Color', plt.color.green, 'LineWidth', 2 )
xlabel('$x$ [m]')
ylabel('energy $e$')
axis( [ st_param.x0, st_param.x5, 1.6, 2.6 ] )
grid on

% % % Save plot to eps file
% fileName = 'plots/shocktube_profiles';
% saveas( gcf, fileName, 'epsc' )
% fprintf('Saved graph to file %s.eps.\n', fileName);

end
