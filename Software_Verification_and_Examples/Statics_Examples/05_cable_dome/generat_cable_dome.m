function [N,C_b,C_s,C] =generat_cable_dome(R,p,m,h,beta)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function gives the configuration of a Levy cable dome
%
% Inputs:
%   R: radius of the cable dome
%   p: complexity in circular direction
%   m: complexity in radial direction
%   h: height of the roof
%   beta: angle of diagonal strings
% Outputs:
%	N: nodal coordinate matrix
%   C_b,C_s: connectivity matrix of bars and strings
%   C: connectivity matrix of all members

%% nodal coordinate
l=2*R;
Rd=(h^2+(l/2)^2)/(2*h); %radius of the arch
%% generate node in one unit
N0=zeros(3,2*m+1);      %initial N
for i=1:m
    N0(:,i)=[i*l/2/(m+1)-R;0;sqrt(Rd^2-((m+1-i)*l/2/(m+1))^2)-(Rd-h)];
end
dy=l/2/(m+1)*tan(beta)+diff([0,N0(3,1:m)])';
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


%% connectivity matrix of a unit
C_b_in_unit = [[1:m]',[m+1:2*m]'];  % Bar 
% Convert above index notation into actual connectivity matrices
C_b2unit = tenseg_ind2C(C_b_in_unit,N);

u=2*m+1;           % # of nodes in a unit
C_s_in_unit1=[1 u   % 2 top 2 diagonal 1 circlur
    1 2*u
    m+1 u
    m+1 2*u
    m+1 m+1+u];%strings of the first loop
C_s_in_temp=[1 2  % strings of inner loop
    1 2+u
    1 m+2
    1 m+2+u
    m+2 m+2+u
    3 2
    3 2+u
    m+3 2
    m+3 2+u
    m+3 m+3+u];%2 top 2 diagonal 1 circlur 2 top 2 diagonal 1 circlur

if rem(m,2)==1     %odd m
C_s_in_unit=[C_s_in_unit1;[kron(2*[0:(m-1)/2-1]',ones(10,2))+kron(ones((m-1)/2,1),C_s_in_temp)]];
C_s_in_unit=[C_s_in_unit;m,m+u];     % top circular strings
else if m==2
    C_s_in_unit=[1 u   % 2 top 2 diagonal 1 circlur
    1 2*u
    m+1 u
    m+1 2*u
    m+1 m+1+u
    1 2  
    1 2+u
    1 m+2
    1 m+2+u
    m+2 m+2+u
    m,m+u];
else if rem(m,2)==0&&m~=2
C_s_in_unit=[C_s_in_unit1;[kron(2*[0:(m-2)/2-2]',ones(10,2)),kron(ones((m-2)/2,1),C_s_in_temp)]];
C_s_in_unit=[C_s_in_unit;           
    m-1 m  
    m-1 m+u
    m-1 2*m
    m-1 2*m+u
    2*m 2*m+u
m,m+u];     % top circular strings
end
    end
end

C_s2unit = tenseg_ind2C(C_s_in_unit,N);
tenseg_plot(N,C_b2unit,C_s2unit)  ;    %plot a unit
title('A unit');
%% replicate connectivity matrix of a unit 
C_in_unit=[C_b_in_unit;C_s_in_unit];
member_unit=size(C_in_unit,1);

C_in=kron(ones(p,1),C_in_unit)+kron([0:p-1]',(2*m+1)*ones(size(C_in_unit))); %C index for the structure

index_n=find(C_in>p*(2*m+1));   %find the index > number of nodes
C_in(index_n)=rem(C_in(index_n),p*(2*m+1));

C = tenseg_ind2C(C_in,N);
ne=size(C,1);%number of element
index_b=kron(ones(p,1),[1:size(C_b_in_unit,1)]')+kron([0:p-1]',member_unit*ones(size(C_b_in_unit,1),1));   
C_b=C(index_b,:);
index_s=setdiff([1:ne],index_b);   
C_s=C(index_s,:);
end

