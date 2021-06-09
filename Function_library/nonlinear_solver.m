function [E,f,q] = nonlinear_solver(data,Xb0,w_t,dXb_t,l0_t,E,p_s,p_e,subsubstep,material)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% This function solves the nonlinear equalibirum equation, considering the
% nonlinear properties of the material. 
% data include the info of structure. 
% p_s, p_e is payload ratio, substep is substep of the payload
% minimize total energy? (1: use, 0: not use) it's time consuming
%%
global Ia Ib C w ne Xb Xa dXa f_int l_int 
use_energy=0;


C=data.C;
ne=data.ne;
Ia=data.Ia;
Ib=data.Ib;
index_b=data.index_b;
index_s=data.index_s;
% E0=data.E;
consti_data=data.consti_data;
A=data.A;


    w_t_m=w_t(:,p_s)+(w_t(:,p_e)-w_t(:,p_s))*linspace(0,1,subsubstep);
    dXb_t_m=dXb_t(:,p_s)+(dXb_t(:,p_e)-dXb_t(:,p_s))*linspace(0,1,subsubstep);    
    l0_t_m=l0_t(:,p_s)+(l0_t(:,p_e)-l0_t(:,p_s))*linspace(0,1,subsubstep);
    
    

% lamda=linspace(p_s,p_e,subsubstep);    %coefficient for substep
for k=1:subsubstep
    w=w_t_m(:,k);               %external force
    Xb=Xb0+dXb_t_m(:,k);         %forced node displacement
    l0=l0_t_m(:,k);            %forced enlongation of string
    disp(['substep',num2str(k)]);
    u=1e-1;
    
        X=[Ia';Ib']\[Xa;Xb];
        l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
        q=E.*A.*(1./l0-1./l);      %force density
       l_int=l;   f_int=q.*l;
    
    for i=1:1e3
        X=[Ia';Ib']\[Xa;Xb];
        l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length
        q=E.*A.*(1./l0-1./l);      %force density
        q_bar=diag(q);
        
        K=kron(C'*q_bar*C,eye(3));                      %stiffness matrix
        Fp=w-K*X;                                       %unbalanced force
        Fp_a=Ia'*Fp;                                 %see the norm of unbalanced force
        if norm(Fp_a)<1e-4
            break
        end        
        N=reshape(X,3,[]);
        B=N*C';
        for j=1:ne
            Ki{j,1}=q_bar(j,j)*eye(3)+E(j)*A(j)*l(j)^(-3)*B(:,j)*B(:,j)';
        end
        K_t=kron(C',eye(3))*blkdiag(Ki{:})*kron(C,eye(3));
        K_taa=Ia'*K_t*Ia;
        
        %modify the stiffness matrix        
        [V_mode,D]=eig(K_taa);                       %刚度矩阵特征根
        d=diag(D);                            %eigen value        
        lmd=min(d);                     %刚度矩阵最小特征根
        if lmd>0
            Km=K_taa+u*eye(size(K_taa)); %修正的刚度矩阵
        else
            Km=K_taa+(abs(lmd)+u)*eye(size(K_taa));
        end
        dXa=Km\Fp_a;
        
        x=1;
%         line search
        if use_energy==1
            opt=optimset('TolX',1e-5);
            [x,V]=fminbnd(@energy,0,1,opt);
        end        
        Xa=Xa+x*dXa;
    end
    
      % change youngs mudulus if string slack
    strain=(l-l0)./l0;        %strain of member    
    [E,sigma]=stress_strain(consti_data,index_b,index_s,strain,material);
    f=sigma.*A;         %member force
    q=f./l;      %reculate force density

end
end
