function tenseg_plot_exforce(Ib,tspan,w_t,node_w,dnb_t,dnb_d_t,dnb_dd_t,node_nb,saveimg)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the acceleration, velocity, displacement of ground
% motion.
%
% Inputs:
%   Ib:transform matrix of pinned coordinates
%   tspan: time sequence
%   w_t: time history of external force
%   node_w: the external force of node_w is ploted
%   dnb_t: displacement of pinned nodes
%   dnb_d_t: velocity of pinned nodes
%   dnb_dd_t: acceleration of pinned nodes
%   saveimg: save image, 1 for yes, 0 for no
% Outputs:
%	a plot of dz_a_t,dz_v_t,dz_d_t
%%
% calculate forced nodal movement of all coordinate:n
dn_t=Ib*dnb_t;
dn_d_t=Ib*dnb_d_t;
dn_dd_t=Ib*dnb_dd_t;

figure
subplot(4,1,1);
plot(tspan,w_t(3*node_w-2,:),tspan,w_t(3*node_w-1,:),tspan,w_t(3*node_w,:),'linewidth',1);
xlabel('t (s)','fontsize',12);
ylabel('F (N)','fontsize',12);
for i = 1:numel(node_w)
    legendInfo{i} = [num2str(node_w(i)) 'X'];
end
for i = 1:numel(node_w)
    legendInfo{end+1} = [num2str(node_w(i)) 'Y'];
end
for i = 1:numel(node_w)
    legendInfo{end+1} = [num2str(node_w(i)) 'Z'];
end
legend(legendInfo,'Interpreter','latex')
legendInfo=[];
subplot(4,1,2);
plot(tspan,dn_t(3*node_nb-2,:),tspan,dn_t(3*node_nb-1,:),tspan,dn_t(3*node_nb,:),'linewidth',1)
xlabel('t (s)','fontsize',12);
ylabel('d (m)','fontsize',12);
for i = 1:numel(node_nb)
    legendInfo{i} = [num2str(node_nb(i)) 'X'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Y'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Z'];
end
legend(legendInfo,'Interpreter','latex')
legendInfo=[];
subplot(4,1,3);
plot(tspan,mean(dn_d_t(1:3:end,:)),tspan,mean(dn_d_t(2:3:end,:)),tspan,mean(dn_d_t(3:3:end,:)),'linewidth',1)
xlabel('t (s)','fontsize',12);
ylabel('v (m/s)','fontsize',12);
for i = 1:numel(node_nb)
    legendInfo{i} = [num2str(node_nb(i)) 'X'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Y'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Z'];
end
legend(legendInfo,'Interpreter','latex')
legendInfo=[];
subplot(4,1,4);
plot(tspan,mean(dn_dd_t(1:3:end,:)),tspan,mean(dn_dd_t(2:3:end,:)),tspan,mean(dn_dd_t(3:3:end,:)),'linewidth',1)
xlabel('t (s)','fontsize',12);
ylabel('a (m/s^2)','fontsize',12);
for i = 1:numel(node_nb)
    legendInfo{i} = [num2str(node_nb(i)) 'X'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Y'];
end
for i = 1:numel(node_nb)
    legendInfo{end+1} = [num2str(node_nb(i)) 'Z'];
end
legend(legendInfo,'Interpreter','latex')
legendInfo=[];
currentFigure = gcf;
title(currentFigure.Children(end), 'Average force and ground motion','fontsize',15);
if saveimg==1
    saveas(gcf,'Average force and ground motion.png');
end
end

%%
% figure
% subplot(4,1,1);
% plot(tspan,mean(w_t(1:3:end,:)),tspan,mean(w_t(2:3:end,:)),tspan,mean(w_t(3:3:end,:)),'linewidth',1)
% xlabel('t (s)','fontsize',18);
% ylabel('F (N)','fontsize',18);
% legend('X','Y','Z');
% subplot(4,1,2);
% plot(tspan,mean(dnb_t(1:3:end,:)),tspan,mean(dnb_t(2:3:end,:)),tspan,mean(dnb_t(3:3:end,:)),'linewidth',1)
% xlabel('t (s)','fontsize',18);
% ylabel('d (m)','fontsize',18);
% legend('X','Y','Z');
% subplot(4,1,3);
% plot(tspan,mean(dnb_d_t(1:3:end,:)),tspan,mean(dnb_d_t(2:3:end,:)),tspan,mean(dnb_d_t(3:3:end,:)),'linewidth',1)
% xlabel('t (s)','fontsize',18);
% ylabel('v (m/s)','fontsize',18);
% legend('X','Y','Z');
% subplot(4,1,4);
% plot(tspan,mean(dnb_dd_t(1:3:end,:)),tspan,mean(dnb_dd_t(2:3:end,:)),tspan,mean(dnb_dd_t(3:3:end,:)),'linewidth',1)
% xlabel('t (s)','fontsize',18);
% ylabel('a (m/s^2)','fontsize',18);
% legend('X','Y','Z');
% currentFigure = gcf;
% title(currentFigure.Children(end), 'Average force and ground motion','fontsize',15);
% if saveimg==1
%     saveas(gcf,'ground_motion_info.png');
% end