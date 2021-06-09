function [A,r]=minimass_buckle_hollow(f,l,t,Eb,c)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% this function is to solve radius of hollow bar with thickness t
%input f force of bar;l length of bar; t thickness of hollow bar
%E_b young's modulus of bar;c is the coefficient of safty
% output A area; r radius
%%
global thick force_b length_bar c_e
thick=t; force_b=f; length_bar=l;c_e=c;
r=fsolve(@area_hollowbar,t*ones(size(l)));    %solve r
A=pi*(r.^2-(r-t).^2);
end

function F=area_hollowbar(r)
% F=pi^3*E_b/4/l^2*(r^4-(r-t)^4)-f;
global thick force_b length_bar Eb c_e
F=pi^3*Eb/4./(length_bar.^2).*(r.^4-(r-thick).^4)-force_b/c_e;
end
