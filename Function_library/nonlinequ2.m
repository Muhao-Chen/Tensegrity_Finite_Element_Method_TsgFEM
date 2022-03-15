function X=nonlinequ2(X0)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%solve nonlinear equilibrium equations using pseudo inverse( tensegrity
%with some constraints but not enough to guarantee positive definite of
%stiffness matrix)
%%
global E_bar A_bar l0_bar C w ne Ia Ib

X=X0;
Xb=Ib'*X;
u=1e-6;
for i=1:1e3
    Xa=Ia'*X;
    l_bar=diag(sqrt(sum((reshape(X,3,[])*C').^2)));  % bar length matrix
    q_bar=E_bar*A_bar*(inv(l0_bar)-inv(l_bar));      % force density
    K=kron(C'*q_bar*C,eye(3));                       % stiffness matrix
    Fp=w-K*X;                                        % unbalanced force
    norm(Ia'*Fp)                                     % see the norm of unbalanced force
    if norm(Ia'*Fp)<u
        break
    end
    %calculate stiffness matrix
    % dqdn=zeros(ne,3*n);
    % for j=1:ne
    %     dqdn(j,:)=X'*kron((C(j,:))'*C(j,:),eye(3));
    % end
    % dqdn=E_bar*A_bar*(l_bar)^-3*dqdn;
    %
    % K_t=kron(C'*q_bar*C,eye(3))+...
    %     kron(C',eye(3))*diag(kron(C,eye(3))*X)*kron(eye(ne),ones(3,1))*dqdn;
    N=reshape(X,3,[]);
    B=N*C';
    for i=1:ne
        Ki{i,1}=q_bar(i,i)*eye(3)+E_bar(i,i)*A_bar(i,i)*l_bar(i,i)^(-3)*B(:,i)*B(:,i)';
    end
    K_t=kron(C',eye(3))*blkdiag(Ki{:})*kron(C,eye(3));

    K_taa=Ia'*K_t*Ia;

    dXa=pinv(K_taa)*(Ia'*Fp);
    Xa=Xa+dXa;
    X=[Ia';Ib']\[Xa;Xb];
end
