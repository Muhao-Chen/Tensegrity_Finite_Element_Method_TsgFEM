function [consti_data,Eb,Es]=material_info()
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
%output information of material: strain-stress, strain-modulus
% bar is Q345 steel, string is steel rope
%% input the info of bars
% strain_b1=[0.001267,0.00157,0.02480,0.033,0.06,0.07,0.08];  % strain of bar Q345 
% stress_b1=1e6*[261,273,288.6,313.6,398.5,0,0];              % stress of bar Q345
strain_b1=[1456e-6,23301e-6];  % strain of bar Q345  双线性
stress_b1=1e6*[300,435];              % stress of bar Q345   双线性

data_b0=[strain_b1;stress_b1];            %material info for strain>0
data_b1=[-fliplr(data_b0),[0;0],data_b0];         %bar info from strain to stress
stress_b=data_b1(2,:);
strain_b=data_b1(1,:);
data_b_E=diff(stress_b)./diff(strain_b);
data_b_strain=strain_b(1:end-1);
data_b2=[data_b_strain;data_b_E];               %bar info from strain to modulus

%% input the info of strings
% strain_s=[-1,0,0.016099,0.02,0.021];   % strain of string 钢丝绳
% stress_s=1e6*[0,0,1223.5,0,0];          % stress of string 钢丝绳
% stress_s=2*1e6*[0,0,790.95,0,0];          % stress of string 钢丝绳
% strain_s=[1,0,0.001267,0.00157,0.02480,0.033,0.06,0.07,0.08];  % strain of bar Q345 
% stress_s=4*1e6*[0,0,261,273,288.6,313.6,398.5,0,0];              % stress of bar Q345
strain_s=16*[-1,0,1456e-6,23301e-6];  % strain of bar Q345  双线性
stress_s=4*1e6*[0,0,300,435];              % stress of bar Q345   双线性

data_s1=[strain_s;stress_s];                        %string info from strain to stress
data_s_E=diff(stress_s)./diff(strain_s);
data_s_strain=strain_s(1:end-1);
data_s2=[data_s_strain;data_s_E];                    %string info from strain to modulus

%% output info
% data_b.data_b1=data_b1;
% data_b.data_b2=data_b2;
% data_s.data_s1=data_s1;
% data_s.data_s2=data_s2;
consti_data.data_b1=data_b1;
consti_data.data_b2=data_b2;
consti_data.data_s1=data_s1;
consti_data.data_s2=data_s2;
%% plot the stress-strain curve
% stress_strain
figure
plot(data_b1(1,:),data_b1(2,:),'k-o','linewidth',1.5)
grid on; 
xlabel('应变','fontsize',14);
ylabel('应力/Pa','fontsize',14);title('杆');
% saveas(gcf,'1压杆本构.png');

figure
plot(data_s1(1,:),data_s1(2,:),'k-o','linewidth',1.5)
grid on;
axis([-0.015,0.025,-max(stress_s),max(stress_s)])
xlabel('应变','fontsize',14);
ylabel('应力/Pa','fontsize',14);
% title('索');
% saveas(gcf,'1拉索本构.png');
%% plot the modulus-strain curve
% E-strain
figure
xi = linspace(2*min(data_b_strain),2*max(data_b_strain),1000);
yi = interp1(data_b_strain, data_b_E, xi,'previous',0);
plot(data_b_strain, data_b_E,'ko', xi, yi,'linewidth',1.5);
xlabel('应变','fontsize',14);
ylabel('弹性模量/Pa','fontsize',14);title('杆');
% saveas(gcf,'1压杆弹性模量.png');

figure
xj = linspace(-1,2*max(data_s_strain),1000);
yj = interp1(data_s_strain, data_s_E, xj,'previous',0);
plot(data_s_strain, data_s_E,'ko', xj, yj,'linewidth',1.5);
% axis([-0.04,0.04,-max(data_s_E),max(data_s_E)])
xlabel('应变','fontsize',14);
ylabel('弹性模量/Pa','fontsize',14);
title('索');
% saveas(gcf,'1拉索弹性模量.png');


Eb=interp1(data_b2(1,:), data_b2(2,:), 0,'previous',0);
Es=interp1(data_s2(1,:), data_s2(2,:), 0,'previous',0);

end

