%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Static calculation of a tensegrity tower %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%%
clc;clear;close all;% global l  Eb Es
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties:'linear_elastic' multielastic plastic
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

substep=20;                                     % load steps
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of double layer prism.
R=50;          %radius
p=12;          %complexity for cable dome
m=2;   %number of circle of the vertical bars
h=0.15*2*R;   %hight of the dome
beta=30*pi/180*ones(m,1);    %all angle of diagonal string
[N,C_b,C_s,C] =generat_cable_dome(R,p,m,h,beta);
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

tenseg_plot(N,C_b,C_s);
title('Cable dome');

%% Boundary constraints
pinned_X=(5:5:60)'; pinned_Y=(5:5:60)'; pinned_Z=(5:5:60)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
gr2={[1];[2];[3 4];[5 6];[7];[8 9];[10 11];[12];[13]}; % group of 1 vertical bar(out) ...2 vertical bar(in) 3 top string 4 diagonal string 5 circlur string 6 top string 7 diagonal string 8 circlur string bottom 9 circluar string top
gr2_t=generate_gr_gr1(gr2,p);
Gp=tenseg_str_gp(gr2_t,C);    %generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=6;                   % number of groups with designed force
fd=20^2/3*1000;                        % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);

% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.2,0.8],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.2,0.8],r_s);
R3Ddata.Nradius=1.2*max(R3Ddata.Bradius)*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Cable dome',R3Ddata);

%% input file of ANSYS
ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),fullfile(savePath,'cable_dome'));

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0; % damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn); % critical damping

%% mode analysis
num_plt=7:9; % mode index to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% external force, forced motion of nodes, shrink of strings
% calculate external force and 
ind_w=[(1:5:56)'*3;(2:5:57)'*3];w=[-4e6*ones(12,1);-2e5*ones(12,1)]; % exert force on top node of cable dome, exertnal node
ind_dnb=[]; dnb0=[];
ind_dl0=[]; dl0=[];   %extent rest length of bar
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);

% case 'force' 'displacement' 'rest_length'

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
% data_out=static_solver2(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;

%% plot member force 
tenseg_plot_result(1:substep,t_t(1:3,:),{'element 1','element 2','element 3'},{'Load step','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*1,3*2],:),{'1Z','2Z'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);

%% Plot final configuration
tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],R3Ddata,l0_t(index_s,end));

%% save output data
if savedata==1
    save (fullfile(savePath,['cable_dome',material{1},'.mat',num2str(material{2})]));
end

%% make video of the dynamic
name=fullfile(savePath,['cable_dome','tf_',num2str(tf),material{1},num2str(material{2})]);
% tenseg_video(n_t,C_b,C_s,[],substep,name,savevideo);
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2});

%% linearized dynaimcs 
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
