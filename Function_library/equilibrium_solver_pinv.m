function data_out=nonlinequ(data)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%solve nonlinear equilibrium equations using modified Newton method
%converge to stable equilibrium
% global E_bar A_bar l0_bar Ia Ib C w ne 
%% input data
C=data.C;
ne=data.ne;
Ia=data.Ia;
Ib=data.Ib;
w=data.w;
E=data.E;
A=data.A;
l0=data.l0;
X0=data.N(:);
%% calculate equilibrium
X=X0;
Xb=Ib'*X;
u=1e-1;
for i=1:1e3
Xa=Ia'*X;

l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length 
q=E.*A.*(1./l0-1./l);      %force density
q_bar=diag(q);

K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix
Fp=w-K*X;                                       %unbalanced force
% norm(Ia'*Fp);                                  %see the norm of unbalanced force
if norm(Ia'*Fp)<1e-6
    break 
end

N=reshape(X,3,[]);
B=N*C';
for i=1:ne
Ki{i,1}=q_bar(i,i)*eye(3)+E(i)*A(i)*l(i)^(-3)*B(:,i)*B(:,i)';
end
K_t=kron(C',eye(3))*blkdiag(Ki{:})*kron(C,eye(3));
K_taa=Ia'*K_t*Ia;
dXa=pinv(K_taa)*(Ia'*Fp);
Xa=Xa+dXa;
X=[Ia';Ib']\[Xa;Xb];

end
%% output data
data_out=data;
data_out.N=reshape(X,3,[]);
data_out.l=l;
data_out.q=q;
data_out.V=energy_cal(data_out);
data_out.Fpn=norm(Ia'*Fp);







 

