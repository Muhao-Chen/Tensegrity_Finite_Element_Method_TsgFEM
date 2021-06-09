function [data_b1,data_b2]=blin_consti_data(Eb,sigma_b)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the bilinear constitutive law given Young's modulus and
% yielding stress, 
%
% Inputs:
%   Eb: tangent Young's modulus in zero strain 
%   sigma_b: yielding stress(the first turing point of stress-strain curve)
% Outputs:
%	data_b1: strain information( tension and compression) 
%	data_b2: stress information( tension and compression) 
%% material
strain_b1=[sigma_b/Eb,2];
stress_b1=[sigma_b,sigma_b];
[data_b1,data_b2,Eb,sigma_b]=point2consti_data(strain_b1,stress_b1);
end

