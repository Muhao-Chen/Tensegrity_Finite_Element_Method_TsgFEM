function [A,r]=minimass_buckle_solid(f,l,Eb,c)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% this function is to solve radius of solid bar 
%input f force of bar;l length of bar; 
%E_b young's modulus of bar;c is the coefficient of safty
% output A area; r radius
%%
A=2*l.*sqrt(f/c/pi/Eb);     
r=sqrt(A/pi);
end

