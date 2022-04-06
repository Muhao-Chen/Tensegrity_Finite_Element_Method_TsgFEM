function [N,C_b,C_s] = three_d_bar(p1,p2,q,angle)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%%
% p1 = [0 0 0];p2 = [4 0 0];q = 1; angle = pi/18;
[N,C_b_n,C_s_n]=threedbar(p1,p2,angle,q);

% save all_test
% [y,N,C_b_in,C_s_in]=convert_index_connectivity([1 0 0 0],pi/6,2);
%
% N = round(N,5);
% C_b_n = round(C_b_n,5);
% C_s_n = round(C_s_n,5);
N = setoff_dup(N,6);
C_b_in = transfer_C_b(N,C_b_n);
C_s_in = transfer_C_b(N,C_s_n);
C_b=tenseg_ind2C(C_b_in,N');
C_s=tenseg_ind2C(C_s_in,N');
N = N';
% tenseg_plot(N',C_b,C_s)
end