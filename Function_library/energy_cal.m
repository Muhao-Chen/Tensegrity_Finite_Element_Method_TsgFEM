function V=energy_cal(data)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% calculate the total energy, x is the cofficient in line search
%% input data
N=data.N;
n=N(:);
C=data.C;
ne=data.ne;
Ia=data.Ia;
Ib=data.Ib;
w=data.w;
E=data.E;
A=data.A;
l0=data.l0;


l=sqrt(sum((N*C').^2))'; %bar length 
V=0.5*(l-l0)'*diag(E.*A./l0)*(l-l0)-w'*n;  %结构总势能
end
