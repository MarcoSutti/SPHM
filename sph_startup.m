% Startup.
% Created:     08.07.2011
% Last change: 08.07.2021

close all; clear; clc;
set( 0, 'defaultAxesTickLabelInterpreter', 'latex' );
set( 0, 'defaultLegendInterpreter',        'latex' );
set( 0, 'defaultTextInterpreter',          'latex' );
set( 0, 'defaultAxesFontSize', 14 );

addpath(genpath('../sph_matlab'))

%--------------------------------------------------------------------------
% Define some more modern colors:
plt.color.yellow = [0.894, 0.671, 0.094];
plt.color.green = [0.667, 0.706, 0.118];
plt.color.gray = [0.325, 0.325, 0.325];
plt.color.darkgray = [ 0.2549, 0.2549, 0.2549 ];
plt.color.blue = [0.239, 0.376, 0.655];
plt.color.red = [0.827, 0.341, 0.306];
plt.color.orange = [0.870, 0.443, 0.137];
plt.color.gray1 = [ 0.808, 0.808, 0.808 ];
plt.color.gray2 = [ 0.553, 0.557, 0.557 ];
plt.color.gray3 = [ 0.255, 0.255, 0.255 ];
plt.color.Stone = [ 0.349, 0.471, 0.557 ];
plt.color.Teal = [ 0.282, 0.667, 0.678 ];
plt.color.Cherry = [ 0.600, 0.059, 0.008 ];
plt.color.Cerulean = [ 0.016, 0.573, 0.761 ];
plt.color.Sky = [ 0.388, 0.773, 0.855 ];
%--------------------------------------------------------------------------