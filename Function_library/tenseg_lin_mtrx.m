function [A_lin,B_lin]=tenseg_lin_mtrx(C,n,E,A,l0,M,D,Ia,A_1a)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function calculate the linearized dynamic modal linearized model x_dot=A_lin*x+B_lin*u,
% in which x=[dna;dna_d]; u=[dfa;dl0]; dfa is the force in free nodal coordinates
%
% Inputs:
%   C:connectivity matrix
%   N:nodal coordinate matrix
%   E: Young's modulus vector
%   A: cross sectional area vector
%   l0: rest length vector
%   M: mass matrix
%   D: damping matrix
%   Ia: transfer matrix of free coordinates
%   A_1a: equilibrium matrix of free coordinates with force density as
%   variable
% Outputs:
%	A_lin: linearized state matrix
%   B_lin: linearized input matrix
%% A and B linearized model
K_T=tenseg_stiff_matx2(C,n,E,A,l0);
K_Taa=Ia'*K_T*Ia;
D_aa=Ia'*D*Ia;
M_aa=Ia'*M*Ia;
K_l0a=-A_1a*diag(E.*A.*l0.^(-2));     % sensitive matrix of rest length
% K_l0=A_1;                           % sensitive matrix of force density
% linearized model x_dot=Ax+B*u, in which x=[dna;dna_d]; u=[dfa;dl0]; dfa is the force in free nodal coordinates
A_lin=[zeros(size(Ia,2)),eye(size(Ia,2));-M_aa\K_Taa,-M_aa\D_aa];
B_lin=[zeros(size(Ia,2),size(Ia,2)+numel(l0));inv(M_aa),-M_aa\K_l0a];
end

