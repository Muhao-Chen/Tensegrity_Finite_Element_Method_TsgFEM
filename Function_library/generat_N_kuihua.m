function N=generat_N_kuihua(alph,h)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% This function generate nodal coordinate of a unit of cable dome in complexity m
%%
global l m R p

Rd=(h^2+(l/2)^2)/(2*h); %radius of the arch
%% generate node in one unit
N0=zeros(3,2*m+1);      %initial N
for i=1:m
    N0(:,i)=[i*l/2/(m+1)-R;0;sqrt(Rd^2-((m+1-i)*l/2/(m+1))^2)-(Rd-h)];
end
dy=l/2/(m+1)*tan(alph)+diff([0,N0(3,1:m)])';
N0(:,m+1:2*m)=N0(:,1:m);
N0(3,m+1:2*m)=N0(3,m+1:2*m)-dy';
N0(:,end)=[-R;0;0];
%% rotate node in a unit
beta1=pi/p;
T1=[cos(beta1) -sin(beta1) 0
    sin(beta1) cos(beta1) 0
    0 0 1];
for i=1:m    %rotate nodes
    if rem(i,2)==1
    N0(:,i)=T1* N0(:,i);
     N0(:,i+m)=T1* N0(:,i+m);
    end
end

%% replicate units
beta2=2*pi/p;
T2=[cos(beta2) -sin(beta2) 0
    sin(beta2) cos(beta2) 0
    0 0 1];
N=N0;
for i=1:p-1
    N=[N,T2^i*N0];
end
