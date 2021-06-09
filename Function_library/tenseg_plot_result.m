function tenseg_plot_result(out_tspan,data,legend1,label,name,saveimg)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the result of dynamic analysis.
%
% Inputs:
%   data: data to be ploted, can be coordinate, velocity, force, each row
%   is a term changing with time
%	legend: lengend(1),(2)...
%   lable:[xlabel,ylabel];
% Outputs:
%%	plot the results
figure
for i=1:size(data,1)
    %      n='rb';
    %      plot(out_tspan,data(i,:),n(1,i),'linewidth',1);hold on
    plot(out_tspan,data(i,:),'linewidth',2);hold on
end
set(gca,'fontsize',18,'linewidth',1.15);
legend(legend1,'location','best');
ylabel(label(2),'fontsize',18);
xlabel(label(1),'fontsize',18);
if saveimg==1
    saveas(gcf,name);
end
end

