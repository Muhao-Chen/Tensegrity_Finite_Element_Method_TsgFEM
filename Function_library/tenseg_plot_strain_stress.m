function tenseg_video_strain_stress(data_out,ele_num,time)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function plot the stress strain curve with a moving mark changing
% with time
%
% Inputs:
%   data_out: bar material
%   ele_num:element number to be highlighted (a vector)
%   time: a vector, containing the time point plotted in stress strain
%   curve
% Outputs:
%	a video
% Example:
%   tenseg_video_strain_stress(data_out,[1,32,38],100);
%%
color=['y','m','c','r','g','b','w','k'];
name=['stress_strain3',data_out.material{1}];
figure;
%     set(gcf,'Position',get(0,'ScreenSize'));  %full screen
    for i = 1:size(data_out.strain_t,2)
if strcmp(data_out.material{1},'linear_elastic')|strcmp(data_out.material{1},'plastic')|strcmp(data_out.material{1},'multielastic')
    plot(reshape(data_out.strain_t(ele_num,:)',[],1),reshape(data_out.stress_t(ele_num,:)',[],1),'o','MarkerSize',5);
    % for j=1:numel(ele_num)
    %         plot(data_out.strain_t(ele_num(j),:),data_out.stress_t(ele_num(j),:),'r-','linewidth',2);hold on;
    % end
else
    %             plot(data_out.consti_data.data_b1(1,:),data_out.consti_data.data_b1(2,:),'-o','linewidth',2);
    xx=linspace(min(data_out.consti_data.data_b1(1,:)),max(data_out.consti_data.data_b1(1,:)),10000);
    yy=interp1(data_out.consti_data.data_b1(1,:),data_out.consti_data.data_b1(2,:),xx);
    plot(xx,yy,'-o','linewidth',2);
end
hold on;
% axis([min(min(data_out.strain_t(ele_num,:))),max(max(data_out.strain_t(ele_num,:))),min(min(data_out.stress_t(ele_num,:))),max(max(data_out.stress_t(ele_num,:)))]);
% plot(data_out.strain_t(ele_num,i),data_out.stress_t(ele_num,i),'rs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5]);
    end
for i=1:numel(time)
    for j=1:numel(ele_num)
         [~,pic(i)]=min(abs(data_out.out_tspan-time(i)));  
        plot(data_out.strain_t(ele_num(j),pic(i)),data_out.stress_t(ele_num(j),pic(i)),'rs','LineWidth',2,'MarkerSize',18,'MarkerEdgeColor',color(j),'MarkerFaceColor',color(j));
                hold on;
    end
    
end
        set(gcf,'color','w');
%                 set(gca,'fontsize',18);
        set(gca,'fontsize', 18,'linewidth',1.15)
        set(gcf,'position',[100,100,400,350])
        grid on;
        set(gca, 'XGrid','on'); % XÖáµÄÍø¸ñ
        ylabel('Stress (Pa)','fontsize',18);
        xlabel('Strain','fontsize',18);
% end

