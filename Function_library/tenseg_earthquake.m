function [w_t,dnb_t,dnb_d_t,dnb_dd_t]=tenseg_earthquake(G,C,mass,b,dz_d_t,dz_v_t,dz_a_t,move_ground)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the time history of external force and ground motion
%
% Inputs:
%   G: gravity force
%	dz_d_t: displacement of ground motion with time
%	dz_v_t: velocity of ground motion with time
%   dz_a_t: acceleration of ground motion with time
%   move_ground: % for earthquake, use pinned nodes motion(1)
%   or add inertia force in free node(0)
%
% Outputs:
%   w_t: time history of external force
%   dnb_t: displacement of pinned nodes
%   dnb_d_t: velocity of pinned nodes
%   dnb_dd_t: acceleration of pinned nodes
%% earthquake
switch move_ground
    case 1      % using real ground motion of pinned nodes in dynamic simulation
        % add interia force
        w_t=0*-0.5*kron(abs(C)'*mass,dz_a_t)+G;    %add acel to node, not ground movement
        % forced movement of pinned nodes in XYZ direction
        dnb_t=kron(ones(numel(b)/3,1),dz_d_t);              %move boundary nodes
        dnb_d_t=kron(ones(numel(b)/3,1),dz_v_t);    %velocity of moved boundary nodes
        dnb_dd_t=kron(ones(numel(b)/3,1),dz_a_t);   %acceleration of moved boundary nodes
    case 0      % pinned nodes fixed, add inertia force to free nodal coordinates
        % add interia force
        w_t=-0.5*kron(abs(C)'*mass,dz_a_t)+G;    %add acel to node, not ground movement
        % forced movement of pinned nodes in XYZ direction
        dnb_t=0*kron(ones(numel(b)/3,1),dz_d_t);              %move boundary nodes
        dnb_d_t=0*kron(ones(numel(b)/3,1),dz_v_t);    %velocity of moved boundary nodes
        dnb_dd_t=0*kron(ones(numel(b)/3,1),dz_a_t);   %acceleration of moved boundary nodes
end
end