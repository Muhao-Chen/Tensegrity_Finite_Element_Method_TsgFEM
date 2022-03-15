function [l,q,E,f,n,n_d]=tenseg_dyn_info(t,Y,data_in)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% Input:
%   t: current time
%   Y0: current X, Xd values:   Y0=[X;Xd];
% Output:
%   [l,q,E,f,n,n_d]   information of the structure
%%
% global l q E f n n_d
C=data_in.C;
Ia=data_in.Ia;
Ib=data_in.Ib;
index_b=data_in.index_b;
index_s=data_in.index_s;
% w0=data.w;
% dXb=data.dXb;
consti_data=data_in.consti_data;
slack=data_in.slack;
plastic=data_in.plastic;

A=data_in.A;
l0=data_in.l0;
n0=data_in.N(:);
% w_t=data_in.w_t;           % external force
dnb_t=data.dnb_t;       % forced node displacement
dnb_d_t=data.dnb_d_t;    %velocity of forced moved node
% dnb_dd_t=data.dnb_dd_t;    %velocity of forced moved node

% M=data_in.M;
% D=data_in.D;
dt=data_in.dt;

nf=numel(Y)/2;
na=Y(1:nf,:);               %free node cooridnate
na_d=Y(nf+1:end,:);         %free node velocity

% tspan=0:dt:tf;
ind=floor(t/dt)+1;
dnb=dnb_t(:,ind);
% w=w_t(:,ind);
nb_d=dnb_d_t(:,ind); %this is the velocity of fixed node
% nb_dd=dnb_dd_t(:,ind); %this is the acceleration of fixed node

nb=Ib'*n0+dnb;
n=Ia*na+Ib*nb;
n_d=Ia*na_d+Ib*nb_d;
l=sqrt(sum((reshape(n,3,[])*C').^2))'; %bar length
strain=(l-l0)./l0;        %strain of member
[E,sigma]=stress_strain(consti_data,index_b,index_s,strain,slack,plastic);
f=sigma.*A;         %member force
q=f./l;      %reculate force density
%     q_bar=diag(q);
%     K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix
%
% na_dd=(Ia'*M*Ia)\(Ia'*(-K*n-D*n_d+w));      %dynamic equation
% Yd=[na_d;na_dd];
end