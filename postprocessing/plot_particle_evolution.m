function [] = plot_particle_evolution( geom, plt, time )

% function [] = plot_particle_evolution( geom, plt, time )
% Purpose: Plot particles evolution.
% Created:     2011
% Last change: 20.06.2021

%==========================================================================
%               Plot of particle position and velocity field
%==========================================================================
% figure('units','normalized','outerposition',[0 0 1 1])
%
% figure('Position',[0 0 1346 653])
%
% figure(1)
% clf
% hFig = gcf;
% set(hFig,'units','normalized','outerposition',[0 0 1 1])
% clf
if geom.dim == 1
    if plt.type == 1
        % Plot the particle position
        plot( geom.x, 0*geom.x, 'o', 'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', plt.color.Sky, 'MarkerSize', 7 );
        
        if plt.velocity
            hold on
            % Plot the velocity field
            quiver( geom.x, 0*geom.x, geom.v, 0*geom.v, ...
                1.5, 'LineWidth', 1.25, 'Color', plt.color.blue )
        end
    end
    axis equal
    axis( [ -2.5, 1, -0.5, 0.5 ] )
    title( [ 'time = ', num2str(time), ' s' ] );
    
    pause(0.01)
    hold off
    
elseif geom.dim == 2
    if plt.type == 1
        % Plot the particle position
        plot( geom.x(1:geom.nrp, 1), geom.x(1:geom.nrp, 2), 'o', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', plt.color.Sky, ...
            'MarkerSize', 7 );
        
        if plt.velocity
            hold on
            % Plot the velocity field
            quiver( geom.x(:,1), geom.x(:,2), geom.v(:,1), geom.v(:,2), ...
                1.5, 'LineWidth', 1.25, 'Color', plt.color.blue )
        end
    elseif plt.type == 2
        scatter( geom.x(:,1), geom.x(:,2), 50, sqrt( geom.v(:,1).^2 + geom.v(:,2).^2 ), 'filled' )
        %colorbar;
    elseif plt.type == 3
        scatter( geom.x(:,1), geom.x(:,2), 50, flp.p(:), 'filled' )
        %colorbar;
    end
    axis equal
    %     axis( [ geom.x_min - geom.wf_left, geom.x_max + geom.wf_right, ...
    %         geom.y_min - geom.wf_down, geom.y_max + geom.wf_up ] )
    hold on
    
    plot( geom.x_bp, geom.y_bp, 's', 'MarkerEdgeColor', plt.color.gray3, ...
        'MarkerFaceColor', plt.color.gray1, 'MarkerSize', 7 );
    
    
    title( [ 'time = ', num2str(time), ' s' ] );
    
    pause(0.01)
    hold off
    
elseif geom.dim == 3
    if plt.type==1
        plot3( geom.x(1:geom.nrp, 1), geom.x(1:geom.nrp, 2), geom.x(1:geom.nrp, 3), ...
            'o', 'MarkerEdgeColor', plt.color.Sky, ...
            'MarkerFaceColor', plt.color.Sky, 'MarkerSize', 5)
        axis equal
        axis( [ geom.x_min - geom.wf_left, geom.x_max + geom.wf_right, ...
            geom.y_min - geom.wf_down, geom.y_max + geom.wf_up, ...
            geom.z_min - geom.wf_down, geom.z_max + geom.wf_up ] )   % axis limits
        hold on
        
        plot3( geom.x_bp, geom.y_bp, geom.z_bp, 's', 'MarkerEdgeColor', ...
            'k', 'MarkerFaceColor', plt.color.red, 'MarkerSize', 5 )
        
        if plt.velocity
            % Plot the velocity field
            quiver3( geom.x(:,1), geom.x(:,2), geom.x(:,3), geom.v(:,1), ...
                geom.v(:,2), geom.v(:,3) )
        end
    elseif plt.type==2
        scatter3( geom.x(:,1), geom.x(:,2), geom.x(:,3), 10, ...
            sqrt(geom.v(:,1).^2+geom.v(:,2).^2+geom.v(:,3).^2),'filled')
        colorbar;
        axis equal
        axis( [ geom.x_min - geom.wf_left, geom.x_max + geom.wf_right, ...
            geom.y_min - geom.wf_down, geom.y_max + geom.wf_up, ...
            geom.z_min - geom.wf_down, geom.z_max + geom.wf_up ] )   % axis limits
        hold on
    end
    
    title( [ 'time = ', num2str(time), ' s' ] );
    
    pause(0.01)
    hold off
    
    %figure('units','normalized','outerposition',[0 0 1 1])
    %f2=getframe(gcf);       %gets the gcf
    %mov=addframe(mov,f2);   %adds frames to the AVI file
    %close all
end
end