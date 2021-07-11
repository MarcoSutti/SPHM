function [ sph ] = viscosity( geom, sph )

% function [ sph ] = viscosity( geom, sph )
% Purpose: Defines the fluid particle sph.viscosity eta [kg/(m*s)]

% Created:     ??.??.2011
% Last change: 09.06.2021

%   Jun 9, 2021:
%       The logical indexing should give the same result as the older loop.

sph.eta = zeros( geom.tnp, 1 );


sph.eta( abs(geom.part_type) == 2 ) = 1.0e-3;


% for i=1:geom.nrp
%     if abs(geom.part_type(i))==1          %If we are dealing with a gas...
%         sph.eta(i) = 0;               %... gases do not have sph.viscosity
%     elseif abs(geom.part_type(i))==2      %If we are dealing with water...
%         if sph.visc                 %... if we consider water sph.viscosity...
%             sph.eta(i) = 1.0e-3;      %... assume the value of 1.0e-3 kg/(m*s)
%         end
%     end
% end

end