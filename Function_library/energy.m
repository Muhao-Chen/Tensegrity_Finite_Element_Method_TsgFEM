function V=energy(x)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% calculate the total energy, x is the cofficient in line search
%%
global   A E l0 Ia Ib C Xa Xb dXa  w 

X=Ia*(Xa+x*dXa)+Ib*Xb;
l=sqrt(sum((reshape(X,3,[])*C').^2))'; %bar length 
V=0.5*(l-l0)'*diag(E.*A./l0)*(l-l0)-w'*X;  %结构总势能
end
