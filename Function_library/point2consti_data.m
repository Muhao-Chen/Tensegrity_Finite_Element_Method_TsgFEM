function [data_b1,data_b2,Eb,sigma_b]=point2consti_data(strain_b1,stress_b1)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the constitutive law given stress and strain point
% information
%
% Inputs:
%   strain_b1: bar strain information(only positive number)
%   stress_b1: bar stress information(only positive number)
% Outputs:
%	data_b1: strain information( tension and compression)
%	data_b2: stress information( tension and compression)
%	Eb: tangent Young's modulus in zero strain
%	sigma_b: yielding stress(the first turing point of stress-strain curve)
%%
data_b0=[strain_b1;stress_b1];                    % material info for strain>0
data_b1=[-fliplr(data_b0),[0;0],data_b0];         % bar info from strain to stress
stress_b=data_b1(2,:);
strain_b=data_b1(1,:);
data_b_E=diff(stress_b)./diff(strain_b);
data_b_strain=strain_b(1:end-1);
data_b2=[data_b_strain;data_b_E];                 % bar info from strain to modulus

Eb=interp1(data_b2(1,:), data_b2(2,:), 0,'previous',0);
sigma_b=stress_b1(1);
end