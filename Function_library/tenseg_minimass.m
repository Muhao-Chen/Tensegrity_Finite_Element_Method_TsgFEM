function [A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function gives the minimal mass design of tensegrity in yielding and
% buckling constraints, hollow bar with constant thick or solid bar are
% considered
%
% Inputs:
%   t: members' force
%	l: members' length
%   Gp: group matrix of members
%	sigmas,sigmab: yielding stress of strings and bars
%	Eb,Es: Young's modulus of bars and strings
%   index_b, index_s: number of bars and strings
%   c_b,c_s: coefficient of safety factor of bars and strings
%   rho_b,rho_s: density of bars and strings
%   thick: thickness of hollow bars
%   hollow: use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
%
% Outputs:
%	A_b,A_s,A_gp,A: cross sectional area of bars, strings, member in group, all member
%	r_b,r_s,r_gp,radius:radius of bars, strings, members in group,all,member
%   E: Young's modulus vector
%   l0: rest length vector
%   rho: density vector
%   mass: mass vector
%%
switch hollow
    case 1
        [A_b,r_b]=minimass_buckle_hollow(abs(t(index_b)),l(index_b),thick,Eb,c_b);% output area, radius
    case 0
        [A_b,r_b]=minimass_buckle_solid(abs(t(index_b)),l(index_b),Eb,c_b);% output area, radius
end

ne=numel(t);        %number of elements

% cross sectional of strings
A_s=t(index_s)/sigmas/c_s;	% area of string
r_s=sqrt(A_s/pi);               % radius of string

% cross sectional of all members
I3=eye(ne);
Ind_b=I3(:,index_b);            % index matrix for bar
Ind_s=I3(:,index_s);            % index matrix for string
A=[Ind_b,Ind_s]*[A_b;A_s];      % cross sectional area
A_gp=pinv(Gp)*A;
radius=[Ind_b,Ind_s]*[r_b;r_s];	%radius
r_gp=pinv(Gp)*radius;
E=[Ind_b,Ind_s]*[Eb*ones(numel(index_b),1);Es*ones(numel(index_s),1)];        %Young's modulus vector

%  members' force & rest length
l0=E.*A.*l./(t+E.*A);
% density vector 
rho=zeros(ne,1);
rho(index_b)=rho_b;
rho(index_s)=rho_s;
% mass matrix
mass=rho.*A.*l0;           % mass vector
end

