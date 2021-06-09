function C_b=transfer_C_b(N,C_b)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
[a,b,c]=cart2pol(N(:,1),N(:,2),N(:,3));
chara = 1*a + 3* b +17*c;
C_b_in = C_b(:,1:3);
[a_in,b_in,c_in]=cart2pol(C_b_in(:,1),C_b_in(:,2),C_b_in(:,3));
chara_in = 1*a_in + 3* b_in +17*c_in;
C_b_out = C_b(:,4:6);
[a_out,b_out,c_out]=cart2pol(C_b_out(:,1),C_b_out(:,2),C_b_out(:,3));
chara_out = 1*a_out + 3* b_out +17*c_out;
C_b_1 = [];
C_b_2 = [];
% save test_C_b
for i=1:length(C_b_in(:,1))
%     temp1 = find(abs((chara-chara_in(i)))<0.001);
    temp1 = find(chara==chara_in(i));
    temp2 = find(chara==chara_out(i));
%     temp2 = find(abs((chara-chara_out(i)))<0.001)
    C_b_1 = [C_b_1;temp1(1)];
    C_b_2 = [C_b_2;temp2(1)];
end
C_b = [C_b_1 C_b_2];
end