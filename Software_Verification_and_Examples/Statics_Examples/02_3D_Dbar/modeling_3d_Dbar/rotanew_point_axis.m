function P_left = rotanew_point_axis(P,A_nor,theta)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
A_nor = A_nor/(sqrt(A_nor*A_nor'));
Ax = A_nor(1);Px = P(1);
Ay = A_nor(2);Py = P(2);
Az = A_nor(3);Pz = P(3);
% theta = pi*120/180;
% theta_2 = theta*2;
axp = [Ay*Pz - Az*Py,Az*Px - Ax*Pz,Ax*Py - Ay*Px];
P_left= P*cos(theta) + (axp)*sin(theta) + A_nor*(A_nor*P')*(1 - cos(theta));
% P_right = vpa(P*cos(theta_2) + (axp)*sin(theta_2) + A*(A*P')*(1 - cos(theta_2)),10);