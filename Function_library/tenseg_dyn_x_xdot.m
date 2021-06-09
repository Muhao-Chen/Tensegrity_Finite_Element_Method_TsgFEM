function Yd=tenseg_dyn_x_xdot(t,Y,data_in)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% This function returns the x_dot and x_double_dot.
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   Yd=[Xd,Xdd]
%%
global l q E f n n_d l0 A stress strain
C=data_in.C;
Ia=data_in.Ia;
Ib=data_in.Ib;
index_b=data_in.index_b;
index_s=data_in.index_s;
% w0=data.w;
% dXb=data.dXb;
consti_data=data_in.consti_data;
% slack=data_in.slack;
% plastic=data_in.plastic;
% multielastic=data_in.multielastic;
material=data_in.material;

A=data_in.A;
l0=data_in.l0;
n0=data_in.N(:);
w_t=data_in.w_t;           % external force
dnb_t=data_in.dnb_t;       % forced node displacement
dnb_d_t=data_in.dnb_d_t;    %velocity of forced moved node
dnb_dd_t=data_in.dnb_dd_t;    %velocity of forced moved node

M=data_in.M;
D=data_in.D;
dt=data_in.dt;

nf=numel(Y)/2;
na=Y(1:nf,:);               %free node cooridnate
na_d=Y(nf+1:end,:);         %free node velocity
%%
% tspan=0:dt:tf;
ind=floor(t/dt)+1;
% dnb=dnb_t(:,ind);
% nb_d=dnb_d_t(:,ind); %this is the velocity of fixed node
% nb_dd=dnb_dd_t(:,ind); %this is the acceleration of fixed node

% Get current pinned nodes displacement 
if ischar(dnb_t)
    run(dnb_t) % if dnb is a string, run that script
elseif size(dnb_t,2)==1
    dnb = dnb_t; % dnb can be constant
else
    dnb = dnb_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes velocity 
if ischar(dnb_d_t)
    run(dnb_d_t) % if dnb is a string, run that script
elseif size(dnb_d_t,2)==1
    nb_d = dnb_d_t; % dnb can be constant
else
    nb_d = dnb_d_t(:,ind); % or dnb can be time-varying
end

% Get current pinned nodes acceleration 
if ischar(dnb_dd_t)
    run(dnb_dd_t) % if dnb is a string, run that script
elseif size(dnb_dd_t,2)==1
    nb_dd = dnb_dd_t; % dnb can be constant
else
    nb_dd = dnb_dd_t(:,ind); % or dnb can be time-varying
end


%%
% w=w_t(:,ind);
nb=Ib'*n0+dnb;
n=Ia*na+Ib*nb;
n_d=Ia*na_d+Ib*nb_d;

% Get current external forces
if ischar(w_t)
    run(w_t) % if w_t is a string, run that script
elseif size(w_t,2)==1
    w = w_t; % w_t can be constant
else
    w = w_t(:,ind); % or W can be time-varying
end

l=sqrt(sum((reshape(n,3,[])*C').^2))'; %bar length
strain=(l-l0)./l0;        %strain of member
[E,stress]=stress_strain(consti_data,index_b,index_s,strain,material);
f=stress.*A;         %member force
q=f./l;      %reculate force density
q_bar=diag(q);
K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix

na_dd=(Ia'*M*Ia)\(Ia'*(w-M*Ib*nb_dd-K*n-D*n_d));      %dynamic equation
Yd=[na_d;na_dd];
