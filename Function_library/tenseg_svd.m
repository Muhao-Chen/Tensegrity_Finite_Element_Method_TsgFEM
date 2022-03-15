function [U1,U2,V1,V2,S]=tenseg_svd(A_1ag)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function do the singular value decomposition of equilibrium matrix
%
% Inputs:
%	A_1a: equilibrium matrix with constraints, force density as variable,
%	this variable can be changed to  A_1ag, A_2a, A_2ag
%
% Outputs:
%	U1: column space of A_1ag
%	U2: numspace space of A_1ag transpose: mechanism mode
%	V1: row space of A_1ag
%	V2: numspace space of A_1ag: prestress mode
%%
[U,S,V] = svd(A_1ag);
r=rank(A_1ag);                       % rank of (A_1g)
U1=U(:,1:r);U2=U(:,r+1:end);        % U1 is C(A_1g); U2 is N(A_1g') mechanism mode
% S1=S(1:r,1:r);                      % S1 is singular value of A_1g
V1=V(:,1:r);V2=V(:,r+1:end);        % V1 is C(A_1g'); V2 is N(A_1g) self stress mode
end

