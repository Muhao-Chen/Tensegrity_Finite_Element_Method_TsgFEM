function [consti_data,Eb,Es,sigma_b,sigma_s,rho_b,rho_s]=material_lib(bar_material,string_material)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function is a material library of bars and strings.
% Bar material: Steel_Q345, Carbon_Rod, Steel, UHMWPE, Aluminum, Wood.
% String material: Steel_string, Steel, UHMWPE, Aluminum, Rubber_band.
% If this library does not have the material one wants to analysis, please
% make edits to this function directly.
%
% Inputs:
%   b_material: bar material name
%   s_material: string material name
% Outputs:
%	consti_data.data_b1: strain of bar
%	consti_data.data_b2: stress of bar
%	consti_data.data_s1: strain of string
%	consti_data.data_s2: stress of string
%   Eb: Young's modulus of bar
%   Es: Young's modulus of string
%   sigma_b: yielding stress of bar
%   sigma_s: yielding stress of string
%   rho_b: density of bar
%   rho_s: density of string
% Example:
%   [consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
%%
global  Eb Es
switch bar_material
    case 'Steel_Q345'
        %         rho_b = 7870;
        rho_b = 7870*50;    % amplified density
        strain_b1=[1456e-6,23301e-6,1];  % strain of bar Q345
        stress_b1=1e6*[300,435,435];              % stress of bar Q345
        [data_b1,data_b2,Eb,sigma_b]=point2consti_data(strain_b1,stress_b1);
    case 'Q345_blin'
        rho_b = 7870;
        strain_b1=[1456e-6,23301e-6];  % strain of bar Q345 bilinear
        stress_b1=1e6*[300,435];              % stress of bar Q345
        [data_b1,data_b2,Eb,sigma_b]=point2consti_data(strain_b1,stress_b1);
    case 'Carbon_Rod'
        % -------------- Steel Properties for Bar -----------
        % https://gwcomposites.com/carbon-rods/
        Eb = 138e09;
        rho_b = 1500;
        sigma_b = 1.72e09;
        [data_b1,data_b2]=blin_consti_data(Eb,sigma_b);
    case 'Steel'
        % -------------- Steel Properties for Bar -----------
        Eb = 200e09;
        rho_b = 8000;
        %         rho_b = 80000;
        sigma_b = 300e06;
        [data_b1,data_b2]=blin_consti_data(Eb,sigma_b);
    case 'UHMWPE'
        %---------------- UHMWPE Properties for Bar ----------
        Eb = 120e09;
        rho_b = 970;
        sigma_b = 2.7e09;
        [data_b1,data_b2]=blin_consti_data(Eb,sigma_b);
    case 'Aluminum'
        %---------------- Aluminum for Bar -------------------
        Eb = 60e09;
        rho_b = 2700;
        sigma_b = 110e06;
        [data_b1,data_b2]=blin_consti_data(Eb,sigma_b);
    case 'Wood'
        %---------------- Wood for Bar -------------------
        Eb = 8100e06;
        %         rho_b = 2700;
        rho_b = 2700;
        sigma_b = 78e06;
        [data_b1,data_b2]=blin_consti_data(Eb,sigma_b);
    otherwise
        disp('Edit the material database')
end
switch string_material
    case 'Steel_string'
        rho_s = 7870;
        strain_s1=[0.016099,3];  % strain of string
        stress_s1=1e6*[1223.5,1223.5];              % stress of string
        [data_s1,data_s2,Es,sigma_s]=point2consti_data(strain_s1,stress_s1);
    case 'Steel'
        % -------------- Steel Properties for Strings -----------
        Es = 200e09;
        rho_s = 8000;
        sigma_s = 300e06;
        [data_s1,data_s2]=blin_consti_data(Es,sigma_s);
    case 'Q345_blin'
        rho_s = 7870;
        strain_s1=[1456e-6,23301e-6];  % strain of bar Q345 bilinear
        stress_s1=1e6*[300,435];              % stress of bar Q345
        [data_s1,data_s2,Es,sigma_s]=point2consti_data(strain_s1,stress_s1);
    case 'UHMWPE'
        %------------------ UHMWPE Properties for Strings ---------
        Es = 120e09;
        rho_s = 970;
        sigma_s = 2.7e09;
        [data_s1,data_s2]=blin_consti_data(Es,sigma_s);
    case 'Aluminum'
        %---------------- Aluminum for Strings -------------------
        Es = 60e09;
        %         rho_s = 2700;
        rho_s = 2700*50;    % amplified density not real
        sigma_s = 110e06;
        [data_s1,data_s2]=blin_consti_data(Es,sigma_s);
    case 'Rubber_band'
        %---------------- Rubber band for Strings -------------------
        Es = 2e06;
        rho_s = 1700;
        sigma_s = 1e07;
        [data_s1,data_s2]=blin_consti_data(Es,sigma_s);
    otherwise
        disp('Edit the material database')
end
% bar is Q345 steel, string is steel rope
%% output info
consti_data.data_b1=data_b1;
consti_data.data_b2=data_b2;
consti_data.data_s1=data_s1;
consti_data.data_s2=data_s2;
end