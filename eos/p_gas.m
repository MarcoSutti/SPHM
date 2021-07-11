function [ flp ] = p_gas( flp, sph )

% function [ flp ] = p_gas( flp, sph )
% Gamma law equation of state. Subroutine to calculate the pressure and
% speed of sound for air (ideal gas).
% Created:     ??.??.2011
% Last change: 18.08.2019

% rho    : density                                              [in];
% u      : specific internal energy per unit mass [m^2/s^2]     [in];
% p      : pressure                                            [out];
% c      : sound velocity                                      [out];

gamma = 1.4;   % gamma, or adiabatic index, is the ratio of specific heats:
               % specific heat at constant volume over specific heat at
               % constant pressure (c_v/c_p)

% Gamma law for ideal gas
flp.p = (gamma-1) * flp.rho.*sph.e;

% Speed of sound in the gas medium
flp.c = sqrt( (gamma-1) * sph.e );

end