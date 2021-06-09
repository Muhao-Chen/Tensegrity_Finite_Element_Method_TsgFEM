%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Static calculation of a tensegrity tower %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% Steps of static calculation:
% 1.Specify material properties
% 2.Manual nodal position, connectivity specification
% 3.Boundary pinned nodes
% 4.Group of elements
% 5.Prestress design
% 6.Cross sectional area design
% 7.Generate input file of ANSYS (optional)
% 8.Mass matrix and damping matrix and Modal analysis (optional)
% 9.External force, forced motion of nodes, shrink of strings
% 10.Equilibrium calculation
% 11. Plot and make video, output data to TECPLOT(optional) 
clc;clearvars;close all;

% Specify material properties
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties:'linear_elastic'£¬ 'multielastic'£¬ 'plastic'
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;    % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

substep=100;                                     %ºÉÔØ×Ó²½
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=1;              % save image or not (1) yes (0)no
savedata=1;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no

%% N C of the structure
% Manually specify node positions of a tensegrity tower.
R=10; h=30; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle
for i=1:p               % nodal coordinate matrix N
    N(:,i)=R*[cos(2*pi*(i-1)/p),sin(2*pi*(i-1)/p),0];
end
for i=p+1:2*p
    N(:,i)=[R*cos(2*pi*(i-1)/p+beta*pi/180),R*sin(2*pi*(i-1)/p+beta*pi/180),h];
end
for i=2*p+1:3*p
    N(:,i)=[R*cos(2*pi*(i-1)/p+2*beta*pi/180),R*sin(2*pi*(i-1)/p+2*beta*pi/180),2*h];
end

% Manually specify connectivity indices.
C_b_in = [1 5;2 6;3 4;5 9;6 7;4 8];   % This is indicating the bar connection
% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);

% Manually specify connectivity indices.
C_s_in = [4 5;5 6;6 4;7 8;8 9;9 7;1 4;2 5;3 6;4 7;5 8;6 9];  % This is indicating the string connection
% Convert the above matrices into full connectivity matrices.
C_s = tenseg_ind2C(C_s_in,N);

C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
title('Double layer prism');
%% Boundary constraints
pinned_X=(1:3)'; pinned_Y=(1:3)'; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
gr={(1:3);(4:6);(7:9);(10:12);(13:15);(16:18)};     % number of elements in one group
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);
%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);
%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;
%prestress design
index_gp=[1,2];                 % number of groups with designed force
fd=-1e5*ones(2,1);              % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.2,0.8],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.2,0.8],r_s);
R3Ddata.Nradius=ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,10);

%% external force, forced motion of nodes, shrink of strings
ind_w=4*3-2;w=9*1000;
ind_dnb=[3*[7:9]']; dnb0=5*ones(3,1);
ind_dl0=1; dl0=-0.3;
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);
% calculate external force and 
% ind_w=[];w=[];
% ind_dnb=1*3-2; dnb0=0.5;
% ind_dl0=[]; dl0=[];
% ind_w=[3*[4:6]'-2];w=[500*ones(3,1)];
% ind_dnb=4*3-2; dnb0=0.5;

%% input file of ANSYS
ansys_input_gp_00(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),ind_w,w,ind_dnb,dnb0,'tower_ansys');
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'tower');

%% equilibrium calculation
% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia_new; data.Ib=Ib_new;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t;% forced movement of pinned nodes
data.l0_t=l0_t;% rest length of members
data.substep=substep;    % substep
% nonlinear analysis
data_out=static_solver2(data);        %solve equilibrium using mNewton method

% data_out=static_solver(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;

%% plot member force 
tenseg_plot_result(1:substep,t_t([1;7;13],:),{'bar','horizontal string','vertical string'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);
%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);
%% Plot final configuration
tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end))
%% save output data
if savedata==1
    save (['tower_static_',material{1},'.mat']);
end
%% make video of the dynamic
name=['tower_',material{1},'_slack_',num2str(material{2})];
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2})
% tenseg_video(n_t,C_b,C_s,[],min(substep,50),name,savevideo,R3Ddata);


%% linearized dynaimcs 
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);