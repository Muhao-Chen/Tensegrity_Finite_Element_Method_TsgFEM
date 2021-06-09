function [V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,n,E,A,l0,M,num_plt,saveimg,ampli)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function solves the free vibration mode and frequency of tensegrity
% in given shape.
%
% Inputs:
%   Ia: Tranform matrix to get free nodal coordinate: na=Ia'*n
%	C_b: Bar connectivity matrix
%	C_s: String connectivity matrix
%	C: Connectivity matrix of the whole structure
%	n: Nodal coordinate vector:n=N(:)
%   E: Tangent Young's modulus vector
%   A: Cross sectional area vector
%   l0: Rest length vector
%   M: Mass matrix
%   num_plt: number of modes to be ploted
%   saveimg: ave image or not (1) yes (0)no
%   ampli: coefficient of amplification for deformed shape from initial shape
%
% Outputs:
%	V_mode: this is the matrix whose ith vector is the ith modal shape
%	omega: this is the frequency vector in Hz
%%
switch nargin
    case 11
        ampli=0.1 %default value of ampli=0.1
end
%% calculate the mode shape and frequency
N=reshape(n,3,[]);      %nodal coordinate matrix
% K_T=tenseg_stiff_matx2(C,n,E,A,l0,S);
K_T=tenseg_stiff_matx2(C,n,E,A,l0);
[V_mode,D1] = eig(Ia'*K_T*Ia,Ia'*M*Ia);         % calculate vibration mode
d1=diag(D1);                                    % eigen value
omega=real(sqrt(d1))/2/pi;                   % frequency in Hz
%% plot the mode
figure
plot(1:size(Ia,2),omega,'k-o','linewidth',1.5);
set(gca,'fontsize',18);
xlabel('Order of Vibration Mode','fontsize',18,'Interpreter','latex');
ylabel('Frequency (Hz)','fontsize',18,'Interpreter','latex');
% grid on;
if saveimg==1
    saveas(gcf,'frequency.png');
end

for i=1:numel(num_plt)
    f1=figure;
    title=({['mode ',num2str(num_plt(i))];['f=',num2str(omega(num_plt(i)),'%.4f'),'Hz']});
    %plot buckling mode
    tenseg_plot(N+ampli*max(l0)*reshape(Ia*V_mode(:,num_plt(i)),3,[]),C_b,C_s,f1,[],[],title);
    tenseg_plot_dash(N,C_b,C_s,f1,[],[],title);
    % axis off;
    %     view(2);
    if saveimg==1
        saveas(gcf,['Mode',num2str(num_plt(i)),'.png']);
    end
end

%%  plot 3 mode in one figure
% f2=figure;
% for i=1:numel(plt)
%     %plot buckling mode
%     subplot(1,numel(plt),i);
%     title=(['f = ',num2str(d_omega1(plt(i)),'%.4f'),' Hz']);
%     %     title=[];
%     tenseg_plot(N-10*max(l)*reshape(Ia*V_mode1(:,plt(i)),3,[]),C_b,C_s,f2,[],[],title);
%     xlim([0 11]); ylim([-0.5 1.5]); axis equal;
%     tenseg_plot_dash(N,C_b,C_s,f2,[],[],title);
%     xlim([0 11]); ylim([-0.5 1.5]); axis equal;
%     set(gca,'fontsize',18);
%     %         xlabel('Mode Order');ylabel('Displacement (m)')
%     axis off
%     %         xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% end
% saveas(gcf,['Mode_tower_first4','.png']);
end

