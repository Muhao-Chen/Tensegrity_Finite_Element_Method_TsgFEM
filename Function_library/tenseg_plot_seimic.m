function tenseg_plot_seimic(tspan,dz_a_t,dz_v_t,dz_d_t,saveimg);
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the acceleration, velocity, displacement of ground
% motion
%
% Inputs:
%   tspan: time sequence
%   dz_a_t: acceleration
%   dz_v_t: velocity
%   dz_d_t: displacement
% Outputs:
%	a plot of dz_a_t,dz_v_t,dz_d_t
%%
figure
subplot(3,1,1);
plot(tspan,dz_a_t,'linewidth',1)
xlabel('t (s)','fontsize',18);
ylabel('a (m/s^2)','fontsize',18);
legend('X','Y','Z');
subplot(3,1,2);
plot(tspan,dz_v_t,'linewidth',1)
xlabel('t (s)','fontsize',18);
ylabel('v (m/s)','fontsize',18);
legend('X','Y','Z');
subplot(3,1,3);
plot(tspan,dz_d_t,'linewidth',1)
xlabel('t (s)','fontsize',18);
ylabel('d (m)','fontsize',18);
legend('X','Y','Z');
if saveimg==1
    saveas(gcf,'ground_motion_info.png');
end
end
