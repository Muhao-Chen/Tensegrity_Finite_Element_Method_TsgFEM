function plot_data_torus(data_out,direction,saveimg)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
%plot the information of torus R=10
%%
N_out=data_out.N_out;
t_step=data_out.t_out;          %member force in every step
n_step=data_out.n_out;          %nodal coordinate in every step
substep=size(n_step,2);
%% plot member force and node coordinate

% plot_data_torus(data_out,direction,saveimg);
zb=1:substep;
if direction==1
    figure
    plot(zb,t_step(147,:),'k-o',zb,t_step(152,:),'k-^','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外脊索','内脊索')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1脊索内力.png');
    end
    
    figure
    plot(zb,t_step(149,:),'k-o',zb,t_step(154,:),'k-^','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外斜索','内斜索','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1斜索内力.png');
    end
    
    figure
    plot(zb,t_step(372,:),'k-o',zb,t_step(156,:),'k-^',zb,t_step(157,:),'k-v','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外环索','内环索','内环顶索','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环索内力.png');
    end
    
    figure
    plot(zb,t_step(145,:),'k-o',zb,t_step(146,:),'k-^','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外竖杆','内竖杆')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1竖杆内力.png');
    end
    
    figure
    plot(zb,t_step(4,:),'k-o',zb,t_step(19,:),'k-^','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('水平杆','稳定杆')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环杆内力.png');
    end
    
    figure
    plot(zb,t_step(55,:),'k-o',zb,t_step(76,:),'k-^','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('环竖索','环斜索')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环的索内力.png');
    end
    
    z_whd_tv=n_step(55*3,:)-n_step(55*3,1);
    z_nhd_tv=n_step(56*3,:)-n_step(56*3,1);
    z_whd_bv=n_step(57*3,:)-n_step(57*3,1);
    z_nhd_bv=n_step(58*3,:)-n_step(58*3,1);
    figure
    plot(zb,z_whd_tv,'k-^',zb,z_nhd_tv,'k-o',zb,z_whd_bv,'k-v',zb,z_nhd_bv,'k-*','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外环顶节点（V）','内环顶节点（V）','外环底节点（V）'...
        ,'内环底节点（V）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1节点位移(V).png');
    end
    
    h_whd_th=sqrt(sum((n_step(55*3-[2,1],:)).^2))-norm(n_step(55*3-[2,1],1));
    h_nhd_th=sqrt(sum((n_step(56*3-[2,1],:)).^2))-norm(n_step(56*3-[2,1],1));
    h_whd_bh=sqrt(sum((n_step(57*3-[2,1],:)).^2))-norm(n_step(57*3-[2,1],1));
    h_nhd_bh=sqrt(sum((n_step(58*3-[2,1],:)).^2))-norm(n_step(58*3-[2,1],1));
    figure
    plot(zb,h_whd_th,'k-^',zb,h_nhd_th,'k-o',zb,h_whd_bh,'k-v',zb,h_nhd_bh,'k-*','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('外环顶节点（H）','内环顶节点（H）','外环底节点（H）'...
        ,'内环底节点（H）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1节点位移(H).png');
    end
    
    wdd_h1=sqrt(sum((n_step(2*3-[2,1],:)).^2))-norm(n_step(2*3-[2,1],1));
    wdd_v1=n_step(2*3,:)-n_step(2*3,1);
    wdd_h2=sqrt(sum((n_step(5*3-[2,1],:)).^2))-norm(n_step(5*3-[2,1],1));
    wdd_v2=n_step(5*3,:)-n_step(5*3,1);
    wdd_h3=sqrt(sum((n_step(8*3-[2,1],:)).^2))-norm(n_step(8*3-[2,1],1));
    wdd_v3=n_step(8*3,:)-n_step(8*3,1);
        figure
    plot(zb,wdd_h1,'k-^',zb,wdd_v1,'k-*','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('环外节点（H）','环外节点（V）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环形张拉整体外节点位移.png');
    end
    
    l_14=sqrt(sum((n_step(28*3-[2,1,0],:)-n_step(1*3-[2,1,0],:)).^2))-2*R;
    l_25=sqrt(sum((n_step(34*3-[2,1,0],:)-n_step(7*3-[2,1,0],:)).^2))-2*R;
    l_36=sqrt(sum((n_step(40*3-[2,1,0],:)-n_step(13*3-[2,1,0],:)).^2))-2*R;
    figure
    plot(zb,l_14,'k-o',zb,l_25,'k-^',zb,l_36,'k-v','linewidth',1.5);
    set(gca,'fontsize',18);
    legend('1-10轴支座间距','3-12轴支座间距','5-14轴支座间距','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1支座位移.png');
    end
elseif direction==2
    figure
    plot(zb,t_step(152,:),'k-+',zb,t_step(153,:),'k-o',zb,t_step(178,:),'k-*',zb,t_step(179,:),'k-d',...
        zb,t_step(204,:),'k-.',zb,t_step(205,:),'k-^',zb,t_step(230,:),'k-x',zb,t_step(231,:),'k-v',...
        zb,t_step(256,:),'k->',zb,t_step(257,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内脊索1','内脊索1.5','内脊索3','内脊索3.5','内脊索5','内脊索5.5',...
        '内脊索7','内脊索7.5','内脊索9','内脊索9.5','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环脊索内力.png']);
    end
 
    figure
    plot(zb,t_step(147,:),'k-+',zb,t_step(148,:),'k-o',zb,t_step(173,:),'k-*',zb,t_step(174,:),'k-d',...
        zb,t_step(199,:),'k-.',zb,t_step(200,:),'k-^',zb,t_step(225,:),'k-x',zb,t_step(226,:),'k-v',...
        zb,t_step(251,:),'k->',zb,t_step(252,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外脊索1','外脊索1.5','外脊索3','外脊索3.5','外脊索5','外脊索5.5',...
        '外脊索7','外脊索7.5','外脊索9','外脊索9.5','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环脊索内力.png']);
    end
    
    figure
    plot(zb,t_step(154,:),'k-+',zb,t_step(155,:),'k-o',zb,t_step(180,:),'k-*',zb,t_step(181,:),'k-d',...
        zb,t_step(206,:),'k-.',zb,t_step(207,:),'k-^',zb,t_step(232,:),'k-x',zb,t_step(233,:),'k-v',...
        zb,t_step(258,:),'k->',zb,t_step(259,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内斜索1','内斜索1.5','内斜索3','内斜索3.5','内斜索5','内斜索5.5',...
        '内斜索7','内斜索7.5','内斜索9','内斜索9.5','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环斜索内力.png']);
    end
    
    figure
    plot(zb,t_step(149,:),'k-+',zb,t_step(150,:),'k-o',zb,t_step(175,:),'k-*',zb,t_step(176,:),'k-d',...
        zb,t_step(201,:),'k-.',zb,t_step(202,:),'k-^',zb,t_step(227,:),'k-x',zb,t_step(228,:),'k-v',...
        zb,t_step(253,:),'k->',zb,t_step(254,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外斜索1','外斜索1.5','外斜索3','外斜索3.5','外斜索5','外斜索5.5',...
        '外斜索7','外斜索7.5','外斜索9','外斜索9.5','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环斜索内力.png']);
    end
    
    figure
plot(zb,t_step(157,:),'k-+',zb,t_step(170,:),'k-o',zb,t_step(183,:),'k-*',zb,t_step(196,:),'k-d',...
        zb,t_step(209,:),'k-.',zb,t_step(222,:),'k-^',zb,t_step(235,:),'k-x',zb,t_step(248,:),'k-v',...
        zb,t_step(261,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶索1','内环顶索2','内环顶索3','内环顶索4','内环顶索5','内环顶索6',...
        '内环顶索7','内环顶索8','内环顶索9','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶索内力.png']);
    end
    
    figure
plot(zb,t_step(156,:),'k-+',zb,t_step(169,:),'k-o',zb,t_step(182,:),'k-*',zb,t_step(195,:),'k-d',...
        zb,t_step(208,:),'k-.',zb,t_step(221,:),'k-^',zb,t_step(234,:),'k-x',zb,t_step(247,:),'k-v',...
        zb,t_step(260,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环索1','内环索2','内环索3','内环索4','内环索5','内环索6',...
        '内环索7','内环索8','内环索9','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环索内力.png']);
    end
    
    figure
plot(zb,t_step(372,:),'k-+',zb,t_step(151,:),'k-o',zb,t_step(164,:),'k-*',zb,t_step(177,:),'k-d',...
        zb,t_step(190,:),'k-.',zb,t_step(203,:),'k-^',zb,t_step(216,:),'k-x',zb,t_step(229,:),'k-v',...
        zb,t_step(242,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环索1','外环索2','外环索3','外环索4','外环索5','外环索6',...
        '外环索7','外环索8','外环索9','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环索内力.png']);
    end
    
    figure
plot(zb,t_step(146,:),'k-+',zb,t_step(159,:),'k-o',zb,t_step(172,:),'k-*',zb,t_step(185,:),'k-d',...
        zb,t_step(198,:),'k-.',zb,t_step(211,:),'k-^',zb,t_step(224,:),'k-x',zb,t_step(237,:),'k-v',...
        zb,t_step(250,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内竖杆1','内竖杆2','内竖杆3','内竖杆4','内竖杆5','内竖杆6',...
        '内竖杆7','内竖杆8','内竖杆9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内竖杆内力.png']);
    end
    
    figure
plot(zb,t_step(145,:),'k-+',zb,t_step(158,:),'k-o',zb,t_step(171,:),'k-*',zb,t_step(184,:),'k-d',...
        zb,t_step(197,:),'k-.',zb,t_step(210,:),'k-^',zb,t_step(223,:),'k-x',zb,t_step(236,:),'k-v',...
        zb,t_step(249,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外竖杆1','外竖杆2','外竖杆3','外竖杆4','外竖杆5','外竖杆6',...
        '外竖杆7','外竖杆8','外竖杆9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外竖杆内力.png']);
    end
    
        figure
plot(zb,n_step(3*56,:)-n_step(3*56,1),'k-+',zb,n_step(3*60,:)-n_step(3*60,1),'k-o',zb,n_step(3*64,:)-n_step(3*64,1),'k-*',zb,n_step(3*68,:)-n_step(3*68,1),'k-d',...
        zb,n_step(3*72,:)-n_step(3*72,1),'k-.',zb,n_step(3*76,:)-n_step(3*76,1),'k-^',zb,n_step(3*80,:)-n_step(3*80,1),'k-x',zb,n_step(3*84,:)-n_step(3*84,1),'k-v',...
        zb,n_step(3*88,:)-n_step(3*88,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶节点1（Z）','内环顶节点2（Z）','内环顶节点3（Z）','内环顶节点4（Z）'...
        ,'内环顶节点5（Z）','内环顶节点6（Z）','内环顶节点7（Z）',...
        '内环顶节点8（Z）','内环顶节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*58,:)-n_step(3*58,1),'k-+',zb,n_step(3*62,:)-n_step(3*62,1),'k-o',zb,n_step(3*66,:)-n_step(3*66,1),'k-*',zb,n_step(3*70,:)-n_step(3*70,1),'k-d',...
        zb,n_step(3*74,:)-n_step(3*74,1),'k-.',zb,n_step(3*78,:)-n_step(3*78,1),'k-^',zb,n_step(3*82,:)-n_step(3*82,1),'k-x',zb,n_step(3*86,:)-n_step(3*86,1),'k-v',...
        zb,n_step(3*90,:)-n_step(3*90,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环底节点1（Z）','内环底节点2（Z）','内环底节点3（Z）','内环底节点4（Z）'...
        ,'内环底节点5（Z）','内环底节点6（Z）','内环底节点7（Z）',...
        '内环底节点8（Z）','内环底节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环底节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*55,:)-n_step(3*55,1),'k-+',zb,n_step(3*59,:)-n_step(3*59,1),'k-o',zb,n_step(3*63,:)-n_step(3*63,1),'k-*',zb,n_step(3*67,:)-n_step(3*67,1),'k-d',...
        zb,n_step(3*71,:)-n_step(3*71,1),'k-.',zb,n_step(3*75,:)-n_step(3*75,1),'k-^',zb,n_step(3*79,:)-n_step(3*79,1),'k-x',zb,n_step(3*83,:)-n_step(3*83,1),'k-v',...
        zb,n_step(3*87,:)-n_step(3*87,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环顶节点1（Z）','外环顶节点2（Z）','外环顶节点3（Z）','外环顶节点4（Z）'...
        ,'外环顶节点5（Z）','外环顶节点6（Z）','外环顶节点7（Z）',...
        '外环顶节点8（Z）','外环顶节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环顶节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*57,:)-n_step(3*57,1),'k-+',zb,n_step(3*61,:)-n_step(3*61,1),'k-o',zb,n_step(3*65,:)-n_step(3*65,1),'k-*',zb,n_step(3*69,:)-n_step(3*69,1),'k-d',...
        zb,n_step(3*73,:)-n_step(3*73,1),'k-.',zb,n_step(3*77,:)-n_step(3*77,1),'k-^',zb,n_step(3*81,:)-n_step(3*81,1),'k-x',zb,n_step(3*85,:)-n_step(3*85,1),'k-v',...
        zb,n_step(3*89,:)-n_step(3*89,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环底节点1（Z）','外环底节点2（Z）','外环底节点3（Z）','外环底节点4（Z）'...
        ,'外环底节点5（Z）','外环底节点6（Z）','外环底节点7（Z）',...
        '外环底节点8（Z）','外环底节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环底节点位移Z.png']);
    end
    
    %%%
            figure
plot(zb,n_step(3*56-2,:)-n_step(3*56-2,1),'k-+',zb,n_step(3*60-2,:)-n_step(3*60-2,1),'k-o',zb,n_step(3*64-2,:)-n_step(3*64-2,1),'k-*',zb,n_step(3*68-2,:)-n_step(3*68-2,1),'k-d',...
        zb,n_step(3*72-2,:)-n_step(3*72-2,1),'k-.',zb,n_step(3*76-2,:)-n_step(3*76-2,1),'k-^',zb,n_step(3*80-2,:)-n_step(3*80-2,1),'k-x',zb,n_step(3*84-2,:)-n_step(3*84-2,1),'k-v',...
        zb,n_step(3*88-2,:)-n_step(3*88-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶节点1（X）','内环顶节点2（X）','内环顶节点3（X）','内环顶节点4（X）'...
        ,'内环顶节点5（X）','内环顶节点6（X）','内环顶节点7（X）',...
        '内环顶节点8（X）','内环顶节点9（X）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*58-2,:)-n_step(3*58-2,1),'k-+',zb,n_step(3*62-2,:)-n_step(3*62-2,1),'k-o',zb,n_step(3*66-2,:)-n_step(3*66-2,1),'k-*',zb,n_step(3*70-2,:)-n_step(3*70-2,1),'k-d',...
        zb,n_step(3*74-2,:)-n_step(3*74-2,1),'k-.',zb,n_step(3*78-2,:)-n_step(3*78-2,1),'k-^',zb,n_step(3*82-2,:)-n_step(3*82-2,1),'k-x',zb,n_step(3*86-2,:)-n_step(3*86-2,1),'k-v',...
        zb,n_step(3*90-2,:)-n_step(3*90-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环底节点1（X）','内环底节点2（X）','内环底节点3（X）','内环底节点4（X）'...
        ,'内环底节点5（X）','内环底节点6（X）','内环底节点7（X）',...
        '内环底节点8（X）','内环底节点9（X）','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环底节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*55-2,:)-n_step(3*55-2,1),'k-+',zb,n_step(3*59-2,:)-n_step(3*59-2,1),'k-o',zb,n_step(3*63-2,:)-n_step(3*63-2,1),'k-*',zb,n_step(3*67-2,:)-n_step(3*67-2,1),'k-d',...
        zb,n_step(3*71-2,:)-n_step(3*71-2,1),'k-.',zb,n_step(3*75-2,:)-n_step(3*75-2,1),'k-^',zb,n_step(3*79-2,:)-n_step(3*79-2,1),'k-x',zb,n_step(3*83-2,:)-n_step(3*83-2,1),'k-v',...
        zb,n_step(3*87-2,:)-n_step(3*87-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环顶节点1（X）','外环顶节点2（X）','外环顶节点3（X）','外环顶节点4（X）'...
        ,'外环顶节点5（X）','外环顶节点6（X）','外环顶节点7（X）',...
        '外环顶节点8（X）','外环顶节点9（X）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环顶节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*57-2,:)-n_step(3*57-2,1),'k-+',zb,n_step(3*61-2,:)-n_step(3*61-2,1),'k-o',zb,n_step(3*65-2,:)-n_step(3*65-2,1),'k-*',zb,n_step(3*69-2,:)-n_step(3*69-2,1),'k-d',...
        zb,n_step(3*73-2,:)-n_step(3*73-2,1),'k-.',zb,n_step(3*77-2,:)-n_step(3*77-2,1),'k-^',zb,n_step(3*81-2,:)-n_step(3*81-2,1),'k-x',zb,n_step(3*85-2,:)-n_step(3*85-2,1),'k-v',...
        zb,n_step(3*89-2,:)-n_step(3*89-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环底节点1（X）','外环底节点2（X）','外环底节点3（X）','外环底节点4（X）'...
        ,'外环底节点5（X）','外环底节点6（X）','外环底节点7（X）',...
        '外环底节点8（X）','外环底节点9（X）','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环底节点位移X.png']);
    end
    %%%
    figure
plot(zb,t_step(55,:),'k-+',zb,t_step(56,:),'k-o',zb,t_step(57,:),'k-*',zb,t_step(58,:),'k-d',...
        zb,t_step(59,:),'k-.',zb,t_step(60,:),'k-^',zb,t_step(61,:),'k-x',zb,t_step(62,:),'k-v',...
        zb,t_step(63,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('环竖索1','环竖索2','环竖索3','环竖索4','环竖索5','环竖索6',...
        '环竖索7','环竖索8','环竖索9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'环竖索内力.png']);
    end
    
    figure
plot(zb,t_step(73,:),'k-+',zb,t_step(77,:),'k-o',zb,t_step(81,:),'k-*',zb,t_step(85,:),'k-d',...
        zb,t_step(89,:),'k-.',zb,t_step(93,:),'k-^',zb,t_step(97,:),'k-x',zb,t_step(101,:),'k-v',...
        zb,t_step(105,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('环斜索1','环斜索2','环斜索3','环斜索4','环斜索5','环斜索6',...
        '环斜索7','环斜索8','环斜索9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'环斜索内力.png']);
    end
    
    figure
plot(zb,t_step(2,:),'k-+',zb,t_step(3,:),'k-o',zb,t_step(4,:),'k-*',zb,t_step(5,:),'k-d',...
        zb,t_step(6,:),'k-.',zb,t_step(7,:),'k-^',zb,t_step(8,:),'k-x',zb,t_step(9,:),'k-v',...
        zb,t_step(10,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('水平杆1','水平杆2','水平杆3','水平杆4','水平杆5','水平杆6',...
        '水平杆7','水平杆8','水平杆9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'水平杆内力.png']);
    end
    
    figure
plot(zb,t_step(20,:),'k-+',zb,t_step(22,:),'k-o',zb,t_step(24,:),'k-*',zb,t_step(26,:),'k-d',...
        zb,t_step(28,:),'k-.',zb,t_step(30,:),'k-^',zb,t_step(32,:),'k-x',zb,t_step(34,:),'k-v',...
        zb,t_step(36,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('稳定杆1','稳定杆2','稳定杆3','稳定杆4','稳定杆5','稳定杆6',...
        '稳定杆7','稳定杆8','稳定杆9','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'稳定杆内力.png']);
    end
    
    l_1=sqrt(sum((n_step(28*3-[2,1,0],:)-n_step(1*3-[2,1,0],:)).^2))-2*R;
    l_3=sqrt(sum((n_step(34*3-[2,1,0],:)-n_step(7*3-[2,1,0],:)).^2))-2*R;
    l_5=sqrt(sum((n_step(40*3-[2,1,0],:)-n_step(13*3-[2,1,0],:)).^2))-2*R;
    l_7=sqrt(sum((n_step(46*3-[2,1,0],:)-n_step(19*3-[2,1,0],:)).^2))-2*R;
    l_9=sqrt(sum((n_step(52*3-[2,1,0],:)-n_step(25*3-[2,1,0],:)).^2))-2*R;
      figure
    plot(zb,l_1,'k-o',zb,l_3,'k-^',zb,l_5,'k-v', zb,l_7,'k-.',zb,l_9,'k-x','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('1-10轴支座间距','3-12轴支座间距','5-14轴支座间距','7-16轴支座间距','9-18轴支座间距','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1支座位移.png');
    end
    
    wdd_h1=sqrt(sum((n_step(2*3-[2,1],:)).^2))-norm(n_step(2*3-[2,1],1));
    wdd_v1=n_step(2*3,:)-n_step(2*3,1);
    wdd_h3=sqrt(sum((n_step(8*3-[2,1],:)).^2))-norm(n_step(8*3-[2,1],1));
    wdd_v3=n_step(8*3,:)-n_step(8*3,1);
    wdd_h5=sqrt(sum((n_step(14*3-[2,1],:)).^2))-norm(n_step(14*3-[2,1],1));
    wdd_v5=n_step(14*3,:)-n_step(14*3,1);
    wdd_h7=sqrt(sum((n_step(20*3-[2,1],:)).^2))-norm(n_step(20*3-[2,1],1));
    wdd_v7=n_step(20*3,:)-n_step(20*3,1);
    wdd_h9=sqrt(sum((n_step(26*3-[2,1],:)).^2))-norm(n_step(26*3-[2,1],1));
    wdd_v9=n_step(26*3,:)-n_step(26*3,1);   
       
    figure
     plot(zb,wdd_h1,'k-^',zb,wdd_h3,'k-o',zb,wdd_h5,'k-v',...
        zb,wdd_h7,'k-*',zb,wdd_h9,'k-+','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('环外节点1（H）','环外节点3（H）','环外节点5（H）',...
        '环外节点7（H）','环外节点9（H）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环形张拉整体外节点位移（H）.png');
    end
    
        figure
     plot(zb,wdd_v1,'k-^',zb,wdd_v3,'k-o',zb,wdd_v5,'k-v',...
        zb,wdd_v7,'k-*',zb,wdd_v9,'k-+','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('环外节点1（V）','环外节点3（V）','环外节点5（V）',...
        '环外节点7（V）','环外节点9（V）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,'1环形张拉整体外节点位移（V）.png');
    end
    
    
    
    
else
       figure
    plot(zb,t_step(152,:),'k-+',zb,t_step(153,:),'k-o',zb,t_step(178,:),'k-*',zb,t_step(179,:),'k-d',...
        zb,t_step(204,:),'k-.',zb,t_step(205,:),'k-^',zb,t_step(230,:),'k-x',zb,t_step(231,:),'k-v',...
        zb,t_step(256,:),'k->',zb,t_step(257,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内脊索1','内脊索1.5','内脊索3','内脊索3.5','内脊索5','内脊索5.5',...
        '内脊索7','内脊索7.5','内脊索9','内脊索9.5','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环脊索内力.png']);
    end
 
    figure
    plot(zb,t_step(147,:),'k-+',zb,t_step(148,:),'k-o',zb,t_step(173,:),'k-*',zb,t_step(174,:),'k-d',...
        zb,t_step(199,:),'k-.',zb,t_step(200,:),'k-^',zb,t_step(225,:),'k-x',zb,t_step(226,:),'k-v',...
        zb,t_step(251,:),'k->',zb,t_step(252,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外脊索1','外脊索1.5','外脊索3','外脊索3.5','外脊索5','外脊索5.5',...
        '外脊索7','外脊索7.5','外脊索9','外脊索9.5','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环脊索内力.png']);
    end
    
    figure
    plot(zb,t_step(154,:),'k-+',zb,t_step(155,:),'k-o',zb,t_step(180,:),'k-*',zb,t_step(181,:),'k-d',...
        zb,t_step(206,:),'k-.',zb,t_step(207,:),'k-^',zb,t_step(232,:),'k-x',zb,t_step(233,:),'k-v',...
        zb,t_step(258,:),'k->',zb,t_step(259,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内斜索1','内斜索1.5','内斜索3','内斜索3.5','内斜索5','内斜索5.5',...
        '内斜索7','内斜索7.5','内斜索9','内斜索9.5','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环斜索内力.png']);
    end
    
    figure
    plot(zb,t_step(149,:),'k-+',zb,t_step(150,:),'k-o',zb,t_step(175,:),'k-*',zb,t_step(176,:),'k-d',...
        zb,t_step(201,:),'k-.',zb,t_step(202,:),'k-^',zb,t_step(227,:),'k-x',zb,t_step(228,:),'k-v',...
        zb,t_step(253,:),'k->',zb,t_step(254,:),'k-<','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外斜索1','外斜索1.5','外斜索3','外斜索3.5','外斜索5','外斜索5.5',...
        '外斜索7','外斜索7.5','外斜索9','外斜索9.5','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环斜索内力.png']);
    end
    
    figure
plot(zb,t_step(157,:),'k-+',zb,t_step(170,:),'k-o',zb,t_step(183,:),'k-*',zb,t_step(196,:),'k-d',...
        zb,t_step(209,:),'k-.',zb,t_step(222,:),'k-^',zb,t_step(235,:),'k-x',zb,t_step(248,:),'k-v',...
        zb,t_step(261,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶索1','内环顶索2','内环顶索3','内环顶索4','内环顶索5','内环顶索6',...
        '内环顶索7','内环顶索8','内环顶索9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶索内力.png']);
    end
    
    figure
plot(zb,t_step(156,:),'k-+',zb,t_step(169,:),'k-o',zb,t_step(182,:),'k-*',zb,t_step(195,:),'k-d',...
        zb,t_step(208,:),'k-.',zb,t_step(221,:),'k-^',zb,t_step(234,:),'k-x',zb,t_step(247,:),'k-v',...
        zb,t_step(260,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环索1','内环索2','内环索3','内环索4','内环索5','内环索6',...
        '内环索7','内环索8','内环索9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环索内力.png']);
    end
    
    figure
plot(zb,t_step(372,:),'k-+',zb,t_step(151,:),'k-o',zb,t_step(164,:),'k-*',zb,t_step(177,:),'k-d',...
        zb,t_step(190,:),'k-.',zb,t_step(203,:),'k-^',zb,t_step(216,:),'k-x',zb,t_step(229,:),'k-v',...
        zb,t_step(242,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环索1','外环索2','外环索3','外环索4','外环索5','外环索6',...
        '外环索7','外环索8','外环索9','location','northeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环索内力.png']);
    end
    
    figure
plot(zb,t_step(146,:),'k-+',zb,t_step(159,:),'k-o',zb,t_step(172,:),'k-*',zb,t_step(185,:),'k-d',...
        zb,t_step(198,:),'k-.',zb,t_step(211,:),'k-^',zb,t_step(224,:),'k-x',zb,t_step(237,:),'k-v',...
        zb,t_step(250,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内竖杆1','内竖杆2','内竖杆3','内竖杆4','内竖杆5','内竖杆6',...
        '内竖杆7','内竖杆8','内竖杆9','location','southeast')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内竖杆内力.png']);
    end
    
    figure
plot(zb,t_step(145,:),'k-+',zb,t_step(158,:),'k-o',zb,t_step(171,:),'k-*',zb,t_step(184,:),'k-d',...
        zb,t_step(197,:),'k-.',zb,t_step(210,:),'k-^',zb,t_step(223,:),'k-x',zb,t_step(236,:),'k-v',...
        zb,t_step(249,:),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外竖杆1','外竖杆2','外竖杆3','外竖杆4','外竖杆5','外竖杆6',...
        '外竖杆7','外竖杆8','外竖杆9','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('内力/N','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外竖杆内力.png']);
    end
    
        figure
plot(zb,n_step(3*56,:)-n_step(3*56,1),'k-+',zb,n_step(3*60,:)-n_step(3*60,1),'k-o',zb,n_step(3*64,:)-n_step(3*64,1),'k-*',zb,n_step(3*68,:)-n_step(3*68,1),'k-d',...
        zb,n_step(3*72,:)-n_step(3*72,1),'k-.',zb,n_step(3*76,:)-n_step(3*76,1),'k-^',zb,n_step(3*80,:)-n_step(3*80,1),'k-x',zb,n_step(3*84,:)-n_step(3*84,1),'k-v',...
        zb,n_step(3*88,:)-n_step(3*88,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶节点1（Z）','内环顶节点2（Z）','内环顶节点3（Z）','内环顶节点4（Z）'...
        ,'内环顶节点5（Z）','内环顶节点6（Z）','内环顶节点7（Z）',...
        '内环顶节点8（Z）','内环顶节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*58,:)-n_step(3*58,1),'k-+',zb,n_step(3*62,:)-n_step(3*62,1),'k-o',zb,n_step(3*66,:)-n_step(3*66,1),'k-*',zb,n_step(3*70,:)-n_step(3*70,1),'k-d',...
        zb,n_step(3*74,:)-n_step(3*74,1),'k-.',zb,n_step(3*78,:)-n_step(3*78,1),'k-^',zb,n_step(3*82,:)-n_step(3*82,1),'k-x',zb,n_step(3*86,:)-n_step(3*86,1),'k-v',...
        zb,n_step(3*90,:)-n_step(3*90,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环底节点1（Z）','内环底节点2（Z）','内环底节点3（Z）','内环底节点4（Z）'...
        ,'内环底节点5（Z）','内环底节点6（Z）','内环底节点7（Z）',...
        '内环底节点8（Z）','内环底节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环底节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*55,:)-n_step(3*55,1),'k-+',zb,n_step(3*59,:)-n_step(3*59,1),'k-o',zb,n_step(3*63,:)-n_step(3*63,1),'k-*',zb,n_step(3*67,:)-n_step(3*67,1),'k-d',...
        zb,n_step(3*71,:)-n_step(3*71,1),'k-.',zb,n_step(3*75,:)-n_step(3*75,1),'k-^',zb,n_step(3*79,:)-n_step(3*79,1),'k-x',zb,n_step(3*83,:)-n_step(3*83,1),'k-v',...
        zb,n_step(3*87,:)-n_step(3*87,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环顶节点1（Z）','外环顶节点2（Z）','外环顶节点3（Z）','外环顶节点4（Z）'...
        ,'外环顶节点5（Z）','外环顶节点6（Z）','外环顶节点7（Z）',...
        '外环顶节点8（Z）','外环顶节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环顶节点位移Z.png']);
    end
    

        figure
plot(zb,n_step(3*57,:)-n_step(3*57,1),'k-+',zb,n_step(3*61,:)-n_step(3*61,1),'k-o',zb,n_step(3*65,:)-n_step(3*65,1),'k-*',zb,n_step(3*69,:)-n_step(3*69,1),'k-d',...
        zb,n_step(3*73,:)-n_step(3*73,1),'k-.',zb,n_step(3*77,:)-n_step(3*77,1),'k-^',zb,n_step(3*81,:)-n_step(3*81,1),'k-x',zb,n_step(3*85,:)-n_step(3*85,1),'k-v',...
        zb,n_step(3*89,:)-n_step(3*89,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环底节点1（Z）','外环底节点2（Z）','外环底节点3（Z）','外环底节点4（Z）'...
        ,'外环底节点5（Z）','外环底节点6（Z）','外环底节点7（Z）',...
        '外环底节点8（Z）','外环底节点9（Z）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环底节点位移Z.png']);
    end
    
    %%%
            figure
plot(zb,n_step(3*56-2,:)-n_step(3*56-2,1),'k-+',zb,n_step(3*60-2,:)-n_step(3*60-2,1),'k-o',zb,n_step(3*64-2,:)-n_step(3*64-2,1),'k-*',zb,n_step(3*68-2,:)-n_step(3*68-2,1),'k-d',...
        zb,n_step(3*72-2,:)-n_step(3*72-2,1),'k-.',zb,n_step(3*76-2,:)-n_step(3*76-2,1),'k-^',zb,n_step(3*80-2,:)-n_step(3*80-2,1),'k-x',zb,n_step(3*84-2,:)-n_step(3*84-2,1),'k-v',...
        zb,n_step(3*88-2,:)-n_step(3*88-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环顶节点1（X）','内环顶节点2（X）','内环顶节点3（X）','内环顶节点4（X）'...
        ,'内环顶节点5（X）','内环顶节点6（X）','内环顶节点7（X）',...
        '内环顶节点8（X）','内环顶节点9（X）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环顶节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*58-2,:)-n_step(3*58-2,1),'k-+',zb,n_step(3*62-2,:)-n_step(3*62-2,1),'k-o',zb,n_step(3*66-2,:)-n_step(3*66-2,1),'k-*',zb,n_step(3*70-2,:)-n_step(3*70-2,1),'k-d',...
        zb,n_step(3*74-2,:)-n_step(3*74-2,1),'k-.',zb,n_step(3*78-2,:)-n_step(3*78-2,1),'k-^',zb,n_step(3*82-2,:)-n_step(3*82-2,1),'k-x',zb,n_step(3*86-2,:)-n_step(3*86-2,1),'k-v',...
        zb,n_step(3*90-2,:)-n_step(3*90-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('内环底节点1（X）','内环底节点2（X）','内环底节点3（X）','内环底节点4（X）'...
        ,'内环底节点5（X）','内环底节点6（X）','内环底节点7（X）',...
        '内环底节点8（X）','内环底节点9（X）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'内环底节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*55-2,:)-n_step(3*55-2,1),'k-+',zb,n_step(3*59-2,:)-n_step(3*59-2,1),'k-o',zb,n_step(3*63-2,:)-n_step(3*63-2,1),'k-*',zb,n_step(3*67-2,:)-n_step(3*67-2,1),'k-d',...
        zb,n_step(3*71-2,:)-n_step(3*71-2,1),'k-.',zb,n_step(3*75-2,:)-n_step(3*75-2,1),'k-^',zb,n_step(3*79-2,:)-n_step(3*79-2,1),'k-x',zb,n_step(3*83-2,:)-n_step(3*83-2,1),'k-v',...
        zb,n_step(3*87-2,:)-n_step(3*87-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环顶节点1（X）','外环顶节点2（X）','外环顶节点3（X）','外环顶节点4（X）'...
        ,'外环顶节点5（X）','外环顶节点6（X）','外环顶节点7（X）',...
        '外环顶节点8（X）','外环顶节点9（X）','location','southwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环顶节点位移X.png']);
    end
    

        figure
plot(zb,n_step(3*57-2,:)-n_step(3*57-2,1),'k-+',zb,n_step(3*61-2,:)-n_step(3*61-2,1),'k-o',zb,n_step(3*65-2,:)-n_step(3*65-2,1),'k-*',zb,n_step(3*69-2,:)-n_step(3*69-2,1),'k-d',...
        zb,n_step(3*73-2,:)-n_step(3*73-2,1),'k-.',zb,n_step(3*77-2,:)-n_step(3*77-2,1),'k-^',zb,n_step(3*81-2,:)-n_step(3*81-2,1),'k-x',zb,n_step(3*85-2,:)-n_step(3*85-2,1),'k-v',...
        zb,n_step(3*89-2,:)-n_step(3*89-2,1),'k->','linewidth',1.5);
    set(gca,'fontsize',12);
    legend('外环底节点1（X）','外环底节点2（X）','外环底节点3（X）','外环底节点4（X）'...
        ,'外环底节点5（X）','外环底节点6（X）','外环底节点7（X）',...
        '外环底节点8（X）','外环底节点9（X）','location','northwest')
    xlabel('荷载子步','fontsize',18);
    ylabel('位移/m','fontsize',18);
    if saveimg==1
        saveas(gcf,[num2str(direction),'外环底节点位移X.png']);
    end
 
end

end

