function [C C_b C_s]=generat_C_kuihua(N,m,p)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% This function generates the topology of Kuihua. 
% input N nodal coordinate; m number of loop; p complexity
% output connectivity matrix of Levy dome C, C_b C_s
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