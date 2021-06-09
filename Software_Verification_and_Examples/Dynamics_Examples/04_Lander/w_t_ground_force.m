% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% Example external force script
%
% This script is run at each time step in the simulation, allowing you to
% specify w in terms of n, Ndot,etc

% To illustrate, this is setting an external force on each node that
%    opposes the velocity of the node

%     mb = 1;
%     ms = 0;
%     W = [zeros(2,size(N,2));[-mb/2*9.8*ones(1,12) -ms*9.8*ones(1,size(N,2)-12)]];
%
%     ground_stiffness = 1e4;
%     ground_damping = 3;
%     for position = 1:size(N,2)
%         W1 = zeros(size(N));
%         if N(3,position) <= 0
%             W1(3,position) = -(ground_stiffness*N(3,position)+ground_damping*Nd(3,position));
%             W = W+W1;
%         end
%     end
G=-M*kron(ones(size(M,1)/3,1),[0;0;9.8]);

ground_stiffness=1e6; 
ground_damping = 1e4;
w1=zeros(size(n));
n_z=n(3:3:end);
nd_z=n_d(3:3:end);
n_z_ngtv=n_z;
nd_z_ngtv=nd_z;
n_z_ngtv(find(n_z>0))=0;
nd_z_ngtv(find(n_z>0))=0;
w1(3:3:end)=-(ground_stiffness*n_z_ngtv+ground_damping*nd_z_ngtv);
w=w1+G;


