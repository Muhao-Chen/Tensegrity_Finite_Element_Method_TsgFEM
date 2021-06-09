function [dz_d_t,dz_v_t,dz_a_t]=ground_motion(amplitude,period,tspan,X,Y,Z)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the time history of acceleration,velocity,displacement
% of ground motion, given the amplitude and period of a sinusoidal wave 
%
% Inputs:
%   amplitude: 
%	period
%   X: group motion in X direction 
%	Y: group motion in Y direction
%	Z: group motion in Z direction
%
% Outputs:
%	dz_d_t: displacement of ground motion with time
%	dz_v_t: velocity of ground motion with time
%   dz_a_t: acceleration of ground motion with time
%%
dz_d_t=[1;0;0]*amplitude/(2*pi/period)^2*sin(2*pi/period*tspan);                % displacement of ground motion (time serises)
dz_v_t=-[1;0;0]*amplitude/(2*pi/period)*cos(2*pi/period*tspan);    % velocity of ground motion (time serises)
dz_a_t=[1;0;0]*amplitude*sin(2*pi/period*tspan);    % acceleration of ground motion (time serises)
end

