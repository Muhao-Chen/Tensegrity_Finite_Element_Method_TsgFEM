function [P,P_left_test,P_right_test]=gene_triangle_with_2pt_any(A,B,size)
% this script is used for generate triangle which is vertical to the line
% between these two points
% Given A in the plane
% A and B row in
%%
nmal_vec = A - B;
unit_nmal_vec = nmal_vec /(sqrt(abs(nmal_vec*nmal_vec')));
syms x y z
equal_1=nmal_vec(1)*(x-A(1))+nmal_vec(2)*(y-A(2))...
    +nmal_vec(3)*(z-A(3));
equal_2 = (x-A(1))^2 + (y-A(2))^2+(z-A(3))^2;
equal_3 = (x-A(1))*(B(2)-A(2))+(y-A(2))*(A(1)-B(1));
[x0,y0,z0]=solve(equal_1==0,equal_2==size^2,equal_3==0,x,y,z);
z_index = find(z0==max(z0));
x0 = x0(z_index );
y0 = y0(z_index );
z0 = z0(z_index );
%%
% if A(1)-B(1) ~= 0 && A(2)-B(2) ~= 0 && A(3)-B(3) ~= 0
%     y0 = A(2)+size/sqrt(((A(1)-B(1))/(A(2)-B(2)))^2+1+(((A(1)-B(1))^2+(A(2)-B(2))^2)/((A(2)-B(2))*(A(3)-B(3))))^2);
%     x0 = A(1) + ((A(1)-B(1))/(A(2)-B(2)))*(y0-A(2));
%     z0= A(3) - ((A(1)-B(1))^2+(A(2)-B(2))^2)/((A(2)-B(2))*(A(3)-B(3)))*(y0-A(2));
% elseif A(1)-B(1) == 0 && A(2)-B(2) ~= 0 && A(3)-B(3) ~= 0
%     x0 = A(1);
%     y0 = A(2) + size/sqrt(1+((A(2)-B(2))/((A(3)-B(3))))^2);
%     z0 = A(3) - (A(2)-B(2))/(A(3)-B(3))*(y0-A(2));
% elseif A(1)-B(1) ~= 0 && A(2)-B(2) == 0 && A(3)-B(3) ~= 0
%     y0 = A(2);
%     x0 = A(1) + size/sqrt(1+((A(1)-B(1))/((A(3)-B(3))))^2);
%     z0 = A(3) + (A(1)-B(1))/(A(3)-B(3))*(x0-A(1));
% elseif A(1)-B(1) ~= 0 && A(2)-B(2) ~= 0 && A(3)-B(3) == 0
%     z0 = A(3) + size;
%     y0 = A(2);
%     x0 = A(1);
% elseif A(1)-B(1) == 0 && A(2)-B(2) == 0 && A(3)-B(3) ~= 0
%     z0 = A(3);
%     y0 = A(2) + size/sqrt(2);
%     x0 = A(1) + size/sqrt(2);
% elseif A(1)-B(1) ~= 0 && A(2)-B(2) == 0 && A(3)-B(3) == 0
%     x0 = A(1);
%     y0 = A(2);
%     z0 = A(3) + size;
% elseif A(1)-B(1) == 0 && A(2)-B(2) ~= 0 && A(3)-B(3) == 0
%     y0 = A(2);
%     x0 = A(1);
%     z0 = A(3) + size;
% elseif A(1)-B(1) == 0 && A(2)-B(2) == 0 && A(3)-B(3) == 0
%     disp('You cannot put two same points to generate a structure.')
% end
%%


% x0 = vpa(x0);
% z0 = vpa(z0);
% x0(find(z0~=max(z0)) )= [];
% z0(find(z0~=max(z0)) ) = [];
% if z0(1)<0
%     P = [x0(2) 0 z0(2)];
% else
%     P = [x0(1) 0 z0(1)];
% end
P = [x0 y0 z0];
P = vpa(P);
% A_nor = unit_nmal_vec;
% Ax = A_nor(1);Px = P(1);
% Ay = A_nor(2);Py = P(2);
% Az = A_nor(3);Pz = P(3);
% theta = pi*120/180;
% theta_2 = theta*2;
% axp = [Ay*Pz - Az*Py,Az*Px - Ax*Pz,Ax*Py - Ay*Px];
% P_left= vpa(P*cos(theta) + (axp)*sin(theta) + A*(A*P')*(1 - cos(theta)),10);
% P_right = vpa(P*cos(theta_2) + (axp)*sin(theta_2) + A*(A*P')*(1 - cos(theta_2)),10);
%%
% P_left_test = rotanew_point_axis(P,A_nor,pi/180*120);
% dbstop if error

P_left_test = trans_space_to_rota(P,A,B,pi/180*120);
% P_left_test
% trans_space_to_rota(A,P1,P2,theta)
P_right_test = trans_space_to_rota(P,A,B,pi/180*240);
%
% points = vpa([A;B;P;P_left_test';P_right_test']);
% points
% save point
% plot3(points(:,1),points(:,2),points(:,3),'r*'),axis equal;
