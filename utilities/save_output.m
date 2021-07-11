function [ ] = save_output( x_hist, v_hist, p_hist, e_hist, rho_hist, t_hist, geom, sph, tip )

% function [ ] = save_output( x_hist, v_hist, p_hist, e_hist, rho_hist, t_hist, geom, sph, tip )
% Saves particle information to mat files.
% Created:     20.06.2021
% Last change: 24.06.2021

%   Jun 20, 2021:
%       Filename changes according to example.
%   Jun 20, 2021:
%       Created.

%----------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%----------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName = ['results/example_', num2str(sph.example), '_max_nts_', ...
    num2str(tip.max_nts), '.mat'];

save( fileName, 'x_hist', 'v_hist', 'p_hist', 'e_hist', 'rho_hist', ...
    't_hist', 'geom' )

fprintf('Saved data to file %s.\n', fileName);

end