function P_left = trans_space_to_rota(A,P1,P2,theta)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%%
T = [1 0 0 -P1(1);0 1 0 -P1(2);0 0 1 -P1(3);0 0 0 1];
T_r = [1 0 0 P1(1);0 1 0 P1(2);0 0 1 P1(3);0 0 0 1];
A_t = T*[A 1]';
P1_t = T*[P1 1]';
P2_t = T*[P2 1]';
A_nor = P2_t (1:3)'- P1_t(1:3)';
A_nor = A_nor/(sqrt(A_nor*A_nor'));
P = A_t(1:3)';
P_left = rotanew_point_axis(P,A_nor,theta);
P_left = T_r*[P_left';1];
P_left = vpa(P_left(1:3));
% P_left
end

