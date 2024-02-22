function [ st_results ] = shocktube_exact( st_param )

% Created:     08.09.2022
% Last change: 08.09.2022


st_results.p1 = st_param.pL;
st_results.p5 = st_param.pR;

st_results.rho1 = st_param.rhoL;
st_results.rho5 = st_param.rhoR;

st_results.u1 = st_param.uL;
st_results.u5 = st_param.uR;

%--------------------------------------------------------------------------
gm1 = st_param.gamma - 1;
gp1 = st_param.gamma + 1;

Gamma = gm1/gp1;
beta = gm1/(2*st_param.gamma);

% Handle of the function
fun = @(x) (st_results.p1^beta - x^beta) * sqrt(((1-Gamma^2)*st_results.p1^(1/st_param.gamma))/(Gamma^2 * st_param.rhoL)) - (x - st_results.p5) * sqrt((1-Gamma)/(st_param.rhoR*(x + Gamma * st_results.p5)));

%--------------------------------------------------------------------------
cs1 = sound_speed(st_param.gamma, st_param.pL, st_param.rhoL);

% Velocity of the propagation of the shock:
cs5 = sound_speed(st_param.gamma, st_param.pR, st_param.rhoR);
%--------------------------------------------------------------------------

% st_results.p3 = bisection(fun, 0.2, 0.4, 1e-13)

st_results.p3 = fzero(fun,1);

%--------------------------------------------------------------------------

st_results.p4 = st_results.p3;

st_results.rho4 = st_results.rho5 * (st_results.p4 + Gamma * st_results.p5)/(st_results.p5 + Gamma * st_results.p4);

st_results.u3 = (st_results.p1^beta - st_results.p3^beta) * sqrt(((1-Gamma^2)*st_results.p1^(1/st_param.gamma))/(Gamma^2 * st_param.rhoL));

% Formula equivalente:
% u_3 = u_5 + (P_3 - P_5)/sqrt(rho_5/2 * (gp1 * P_3 + gm1 * P_5))
% u_4 = (P_3 - P_5) * sqrt((1-Gamma)/(rho_R*(P_3 + Gamma * P_5)))

st_results.u4 = st_results.u3;

st_results.rho3 = st_results.rho1 * (st_results.p3/st_results.p1)^(1/st_param.gamma);
%--------------------------------------------------------------------------

cs3 = sound_speed(st_param.gamma, st_results.p3, st_results.rho3);

%--------------------------------------------------------------------------
% Calculate positions:
% Position of the beginning of the rarefaction wave:
st_results.x1 =  st_param.x_barrier - cs1 * st_param.t;                          % 0.26335680867601535;

% Foot of rarefaction
st_results.x2 = st_param.x_barrier + (st_results.u3 - cs3) * st_param.t;                     % 0.4859454374877634;  

% Contact discontinuity:
st_results.x3 = st_param.x_barrier + st_results.u3 * st_param.t;                             % 0.6854905240097902;

% Shock:
st_results.x4 = st_param.x_barrier + cs5 * sqrt(1 + .5 * gp1/st_param.gamma * (st_results.p4/st_results.p5 - 1)) * st_param.t; % 0.8504311464060357;
%--------------------------------------------------------------------------
% Discretization of the intervals:
st_results.x_vec_1 = st_param.x0:(st_results.x1-st_param.x0)/st_param.N:st_results.x1;
st_results.x_vec_2 = st_results.x1:(st_results.x2-st_results.x1)/st_param.N:st_results.x2;
st_results.x_vec_3 = st_results.x2:(st_results.x3-st_results.x2)/st_param.N:st_results.x3;
st_results.x_vec_4 = st_results.x3:(st_results.x4-st_results.x3)/st_param.N:st_results.x4;
st_results.x_vec_5 = st_results.x4:(st_param.x5-st_results.x4)/st_param.N:st_param.x5;
%--------------------------------------------------------------------------

% https://physics.stackexchange.com/questions/423758/how-to-get-exact-solution-to-sod-shock-tube-test
st_results.u2 = 2/gp1 * (cs1 + (st_results.x_vec_2 - st_param.x_barrier)./st_param.t);

factor = 1 - gm1/2 * st_results.u2/cs1;

st_results.rho2 = st_results.rho1 * factor.^(2/gm1);

st_results.p2 = st_results.p1 * factor.^(1/beta);

% Calculate the specific internal energy:
st_results.e1 = st_results.p1./( st_results.rho1 * gm1 );
st_results.e2 = st_results.p2./( st_results.rho2 * gm1 );
st_results.e3 = st_results.p3./( st_results.rho3 * gm1 );
st_results.e4 = st_results.p4./( st_results.rho4 * gm1 );
st_results.e5 = st_results.p5./( st_results.rho5 * gm1 );
% e1 = p1./gm1 + .5 * rho1 .* u1.^2;
% e2 = p2./gm1 + .5 * rho2 .* u2.^2;
% e3 = p3./gm1 + .5 * rho3 .* u3.^2;
% e4 = p4./gm1 + .5 * rho4 .* u4.^2;
% e5 = p5./gm1 + .5 * rho5 .* u5.^2;


fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Positions                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fprintf('Head of rarefaction wave: % .10f\n', st_results.x1 );
fprintf('Foot of rarefaction wave: % .10f\n', st_results.x2 );
fprintf('Contact discontinuity:    % .10f\n', st_results.x3 );
fprintf('Shock:                    % .10f\n', st_results.x4 );

end

%--------------------------------------------------------------------------
function cs = sound_speed(gamma, p, rho)
% Calculates the sound speed.

cs = sqrt(gamma * p/rho);

end
