function K_t=tenseg_stiff_matx2(C,n,E,A,l0,S)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function is to calculate stiffness matrix of tensegrity or clustered
% tensegrity.
%
% Inputs:
%	C: the connectivity matrix
%   n: nodal coordinate vector
%   E: Young's modulus vector
%   A: Cross section vector
%   l0: rest length vector
% Outputs:
%	K_t: the tangent stiffness matrix of the structure
%%
switch nargin
    case 5
        S=eye(size(C,1));
end
%%
N=reshape(n,3,[]);
B=N*C';
l=sqrt(sum(B.^2))';
q=E.*A.*(1./l0-1./l);      %force density
ne=size(C,1);
A1=kron(C',eye(3))*diag(kron(C,eye(3))*n)*...
    kron(eye(ne),ones(3,1));  %equilibrium matrix
Ke=A1*diag(E.*A./(l.^3))*A1';
Kg=kron(C'*diag(q)*C,eye(3));
% K_t=Kg+Ke;         % calculate this way may have round off error:K_t not
% symmetric
K_t=Kg+(Ke+Ke')/2;     % this is to guarantee symmetric real matrix
end