function K_t=tenseg_stiff_matx(data)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%this function is to calculate stiffness matrix of tensegrity
%%
E_bar=diag(data.E);
A_bar=diag(data.A);
l0_bar=diag(data.l0);
C=data.C;
N=data.N;
ne=data.ne;

X=N(:);
l_bar=diag(sqrt(sum((reshape(X,3,[])*C').^2))); %bar length matrix
q_bar=E_bar*A_bar*(inv(l0_bar)-inv(l_bar));      %force density
%calculate stiffness matrix
% dqdn=zeros(ne,3*n);
% for j=1:ne
%     dqdn(j,:)=X'*kron((C(j,:))'*C(j,:),eye(3));
% end
% dqdn=E_bar*A_bar*(l_bar)^-3*dqdn;
%
% K_t=kron(C'*q_bar*C,eye(3))+...
%     kron(C',eye(3))*diag(kron(C,eye(3))*X)*kron(eye(ne),ones(3,1))*dqdn;

% K_taa=Ia'*K_t*Ia;
%%
B=N*C';
K=cell(ne,1);
for i=1:ne
    K{i,1}=q_bar(i,i)*eye(3)+E_bar(i,i)*A_bar(i,i)*l_bar(i,i)^(-3)*B(:,i)*B(:,i)';
end
K_t=kron(C',eye(3))*blkdiag(K{:})*kron(C,eye(3));
end