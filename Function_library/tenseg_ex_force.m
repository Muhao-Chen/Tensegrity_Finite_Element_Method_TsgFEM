function [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,b,type,gravity,acc,C,mass,c_index,amplitude,period)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the time history of external force and ground
% motion.
%
% Inputs:
%   gravity: 1 for considering gravity, 0 for not considering
%   acc: acceleration vector of gravity, for example [0;0;9.8]
%	type: a string containing the type of external force, including:
%           'impulse':exerted to free coordinate
%           'step',
%           'ramp',
%           'vib_force': sinusoidal wave vibration exerted to free coordinate
%           'vib_nodes': sinusoidal wave vibration exerted to pinned coordinate
%   b: index of pinned coordinates
%   c_index: coordinate index in external force or movement, a vector
%   amplitude: amplitude for all type of force, can be a scaler or vector
%	period: only used for sinusoidal wave vibration
%
% Outputs:
%   w_t: time history of external force
%   dnb_t: displacement of pinned nodes
%   dnb_d_t: velocity of pinned nodes
%   dnb_dd_t: acceleration of pinned nodes
%%
switch nargin
    case 9
        period = [];
end

%gravity force vector
G=(gravity)*-0.5*kron(abs(C)'*mass,acc);
%initialize force and displacement data
w_t=zeros(size(G,1),size(tspan,2));     % initialize w_t
dnb_t=zeros(numel(b),numel(tspan));              % move boundary nodes
dnb_d_t=zeros(numel(b),numel(tspan));
dnb_dd_t=zeros(numel(b),numel(tspan));
dz_a_t=[];
% free nodes
% dna_t=zeros(numel(a),numel(tspan));              % free nodes
% dna_d_t=zeros(numel(a),numel(tspan));
% dna_dd_t=zeros(numel(a),numel(tspan));
% n_n = [a;b];
switch type
    case 'impulse'
        w_t(c_index,tspan<0.05)=amplitude*20;        % impluse load in c_index
        w_t=w_t+G*ones(size(tspan));            % add gravity force
    case 'step'
        w_t(c_index,:)=amplitude;        % load in c_index
        w_t=w_t+G*ones(size(tspan));            % add gravity force
    case 'ramp'
        w_t(c_index,:)=amplitude*ones(numel(c_index),1)*linspace(0,1,size(tspan,2));  % load in c_index
        w_t=w_t+G*ones(size(tspan));            % add gravity force
    case 'vib_force'
        %         dz_d_t=-amplitude/(2*pi/period)^2*sin(2*pi/period*tspan);    % displacement of ground motion (time serises)
        %         dz_v_t=-amplitude/(2*pi/period)*cos(2*pi/period*tspan);    % velocity of ground motion (time serises)
        dz_a_t=amplitude*sin(2*pi/period*tspan);    % acceleration of ground motion (time serises)

        w_0=-0.5*kron(abs(C)'*mass,[1;1;1]*dz_a_t);        %load in c_index
        w_t(c_index,:)=w_0(c_index,:);
        w_t=w_t+G*ones(size(tspan));            % add gravity force

        %    [dz_d_t,dz_v_t,dz_a_t]=ground_motion(amplitude,period,tspan,1,0,0);
        %    [w_t,dnb_t,dnb_d_t,dnb_dd_t]=tenseg_earthquake(G,C,mass,b,dz_d_t,dz_v_t,dz_a_t,move_ground);

    case 'vib_nodes'
        w_t=w_t+G*ones(size(tspan));            % add gravity force
        % displacement, velocity, acceleration of ground
        dz_d_t=-amplitude/(2*pi/period)^2*sin(2*pi/period*tspan);    % displacement of ground motion (time serises)
        dz_v_t=-amplitude/(2*pi/period)*cos(2*pi/period*tspan);    % velocity of ground motion (time serises)
        dz_a_t=amplitude*sin(2*pi/period*tspan);    % acceleration of ground motion (time serises)

        dnb_t0=kron(ones(numel(b),1),[1;1;1]*dz_d_t);              %move boundary nodes
        dnb_d_t0=kron(ones(numel(b),1),[1;1;1]*dz_v_t);    %velocity of moved boundary nodes
        dnb_dd_t0=kron(ones(numel(b),1),[1;1;1]*dz_a_t);   %acceleration of moved boundary nodes

        dnb_t(c_index,:)=dnb_t0(c_index,:);
        dnb_d_t(c_index,:)=dnb_d_t0(c_index,:);
        dnb_dd_t(c_index,:)=dnb_dd_t0(c_index,:);

    case 'self_define'
end
end