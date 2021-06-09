%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Jasen linkage with 3 parts%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clc;clearvars;close all;
% global l  Eb Es

%EXAMPLE
clc; clear all; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Wood','Steel_string');
material{1}='linear_elastic'; % index for material properties:'linear_elastic' multielastic plastic
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

substep=100;                                     %荷载子步
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no

%% N C of the structure
% Manually specify node positions (accurate).
N0=[   38.0000   38.0000   -8.7357  -39.6678         0  -19.4476   17.0047   30.3109
    7.8000   22.8000   40.5702   -5.8717         0  -39.6874  -35.4306  -82.5894
         0         0         0         0         0         0         0         0];
     N0(1,:)=N0(1,:)-38;
     N1=N0;N1(3,:)=-5*ones(1,8);
      N2=N0;N2(3,:)=-10*ones(1,8);
      N=[N0,N1,N2];
% N=kron([1 1 1],N0);

% Manually specify connectivity indices.
C_s_in=[];  % This is indicating that string connection
C_b_in0 = [1 2;2 3;3 4;3 5;4 5;4 6;5 7;6 7;6 8;7 8;2 7];  % Similarly, this is saying bar 1 connects node 1 to node 2,
C_b_in=[C_b_in0; C_b_in0+8; C_b_in0+16];   
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=([1 5 9 13 17 21])'; pinned_Y=([1 5 9 13 17 21])'; pinned_Z=(1:nn)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];
Gp=tenseg_str_gp(gr,C);    %generate group matrix 

%% self-stress design 
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=[1];                 % number of groups with designed force
fd=[0];              % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=1:ne;              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
% [A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
A_gp=0.001*ones(ne,1); A=A_gp;
r_b=0.1*ones(ne,1);r_s=[]; r_gp=r_b; radius=r_b;
E=Eb*ones(ne,1);
l0=kron(ones(3,1),[15;50;55.8;41.5;40.1;39.4;39.3;36.7;65.7;49;61.9]);  % rest length should be accurate
% l0=l;
rho=rho_b*ones(ne,1);
mass=A.*l0.*rho;

% Plot the structure with radius
R3Ddata.Bradius=8*r_b;
% R3Ddata.Sradius=r_s;
R3Ddata.Nradius=0.2*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Jasen mechanism',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[];w=[];
ind_dnb=[2*3-2;2*3-1;10*3-2;10*3-1;18*3-2;18*3-1]; dnb0=zeros(6,1);
ind_dl0=[]; dl0=[];
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);
theta1=pi/2; theta2=7/6*pi; theta3=11/6*pi;
dnb_t([find(Ib_new(ind_dnb(1),:)),find(Ib_new(ind_dnb(2),:))],:)=l0(1)*[cos(linspace(0,-2*pi,substep)+theta1)-cos(theta1);sin(linspace(0,-2*pi,substep)+theta1)-sin(theta1)];
dnb_t([find(Ib_new(ind_dnb(3),:)),find(Ib_new(ind_dnb(4),:))],:)=l0(1)*[cos(linspace(0,-2*pi,substep)+theta2)-cos(theta1);sin(linspace(0,-2*pi,substep)+theta2)-sin(theta1)];
dnb_t([find(Ib_new(ind_dnb(5),:)),find(Ib_new(ind_dnb(6),:))],:)=l0(1)*[cos(linspace(0,-2*pi,substep)+theta3)-cos(theta1);sin(linspace(0,-2*pi,substep)+theta3)-sin(theta1)];

%% equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% forced movement of pinned nodes
data.substep=substep;    % substep

% nonlinear analysis
data_out=static_solver(data);        %solve equilibrium using mNewton method

%% Plot final configuration
tenseg_plot_catenary(data_out.N_out{end},C_b,C_s,[],[],[0,90],[],[],l0_t(index_s,end))

t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;

%% input file of ANSYS
ansys_input_gp(N_out{1},C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'Jasen linkage full size');

%% plot member force 
tenseg_plot_result(1:substep,t_t(1:3,:),{'element 1','element 2','element 3'},{'Load step','Force (N)'},'member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*8-2,3*8-1],:),{'8X','8Y'},{'Time (s)','Coordinate (m)'},'X_coordinate.png',saveimg);

%% Plot final configuration
tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end))

%% make video of the dynamic
name=['Jasen_machanism3','tf_',num2str(tf),material{1}];
% tenseg_video(n_t,C_b,C_s,[],substep,name,savevideo);
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,R3Ddata,[0,60],[-80,80,-85,50,-30,30],min(substep,50),name,savevideo,material{2})

%% linearized dynaimcs 
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);

return
%% following code is not used
%% plot member force and node coordinate
figure
plot(1:substep,t_t(1,:),'k-o',1:substep,t_t(2,:),'k-^',1:substep,t_t(3,:),'k-v','linewidth',1.5);
legend('压杆','1索','2索')
xlabel('荷载子步','fontsize',14);
ylabel('内力/N','fontsize',14)
% saveas(gcf,'1环索内力.png');

z_whd1=n_t(10,:)-n_t(10,1);
z_nhd1=n_t(11,:)-n_t(11,1);
figure
plot(1:substep,z_whd1,'k-o',1:substep,z_nhd1,'k-^','linewidth',1.5);
legend('4X','4Y')
xlabel('荷载子步','fontsize',14);
ylabel('位移/m','fontsize',14)
% saveas(gcf,'1节点位移.png');


%% plot structure configuration
tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_t(index_s,end))
% saveas(gcf,'1最终形态.png');

figure(99);
name=['half_Tbar_'];
set(gcf,'Position',get(0,'ScreenSize'));
for n = 1:substep
    tenseg_plot_catenary( N_out{n},C_b,C_s,99,[],[0,0],[],[],l0_t(index_s,n));hold on
    xlim([-1,1]); ylim([-3.5,3.5]); zlim([-1,1]);
    tenseg_savegif_forever(name);
    hold off;
end
close 

