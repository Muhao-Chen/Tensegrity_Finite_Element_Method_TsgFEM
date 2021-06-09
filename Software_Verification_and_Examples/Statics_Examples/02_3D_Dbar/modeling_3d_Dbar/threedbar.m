function [N,C_b_n,C_s_n] = threedbar(p1,p2,angle,q,N,C_b_n,C_s_n)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% p is 1*3
switch nargin
    case 6
        C_s_n = [];
    case 5
        C_s_n  = [];
        C_b_n  = [];
	case 4
		C_s_n  = [];
        C_b_n  = [];
		N = [];
end
mid = (p1+p2)/2;
size_r = tan(angle) * eulerdst(p1,p2)/2;
[up_p,left_p,right_p]=gene_triangle_with_2pt_any(mid,p2,size_r);
up_p = double(up_p);
left_p = double(left_p);
right_p = double(right_p);
% up_p is 1*3,left_p is 3*1,right_p is 3*1
left_p = left_p';
right_p = right_p';
N = [N;p1;p2;up_p;left_p;right_p];
[a,b] = size(C_b_n);
if (a+b)~=0
index_contain_p1_p2 = ismember(C_b_n,[p1 p2],'rows');
C_b_n(index_contain_p1_p2,:) = [];
end
% C_s_n = [C_b_n;p1 up_p;p1 right_p;p1 left_p;p2 up_p;p2 right_p;p2 left_p;p1 mid;mid p2;up_p left_p;left_p right_p;right_p up_p];
% C_b_n = [C_s_n;up_p mid;right_p mid;left_p mid];
C_b_n = [C_b_n;p1 up_p;p1 right_p;p1 left_p;p2 up_p;p2 right_p;p2 left_p];
C_s_n = [C_s_n;up_p left_p;left_p right_p;right_p up_p;p1 p2];
if q==1
    return
else
    [N,C_b_n,C_s_n]=threedbar(p1,up_p,angle,q-1,N,C_b_n,C_s_n);
    [N,C_b_n,C_s_n]=threedbar(p2,up_p,angle,q-1,N,C_b_n,C_s_n);
    [N,C_b_n,C_s_n]=threedbar(p1,right_p,angle,q-1,N,C_b_n,C_s_n);
    [N,C_b_n,C_s_n]=threedbar(p2,right_p,angle,q-1,N,C_b_n,C_s_n);
    [N,C_b_n,C_s_n]=threedbar(p1,left_p,angle,q-1,N,C_b_n,C_s_n);
    [N,C_b_n,C_s_n]=threedbar(p2,left_p,angle,q-1,N,C_b_n,C_s_n);
end
