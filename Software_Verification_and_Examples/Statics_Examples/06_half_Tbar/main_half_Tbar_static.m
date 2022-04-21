%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Self eqilibrium cable dome with Tensegrity torus%%%%%%%%%%%%%
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
%%
clc;clear;close all;
% Specify material properties
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Q345_blin','Q345_blin');
material{1}='linear_elastic'; % index for material properties:linear_elastic; multielastic; plastic
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
thick=6e-3;        % thickness of hollow bar
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

substep=100;                                     % load steps
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions (in 3D).
n1 = [-1 0 0]';n2 = [0 0 0]';
n3 = [1 0 0]'; n4 = [0 0 1]';

% Put node vectors in node matrix. Node matrix has to be 3xn for n nodes.
N = [n1 n2 n3 n4 ];
% Manually specify connectivity indices.
C_s_in=[1 4;3 4];  % This is indicating that string connection
C_b_in = [2 4];  % Similarly, this is saying bar 1 connects node 1 to node 2,

% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);% ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
% tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=(1:3)'; pinned_Y=(1:4)'; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];                     %group index, it is a cell containing vectors, each vector contains # of members in the same group
Gp=tenseg_str_gp(gr,C);    %generate group matrix
n_g=size(Gp,2);            %number of group for elements

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
index_gp=1;                 % number of groups with designed force
fd=-50;              % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional area design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.002,0.008],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.002,0.008],r_s);
R3Ddata.Nradius=0.01*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'2 string 1 bar',R3Ddata);

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% modal analysis
num_plt=1:2;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,0.02);

%% external force, forced motion of nodes, shrink of strings
ind_w=4*3-2;w=10000;
% ind_dnb=4*3-2; dnb0=0.3;
ind_dnb=[]; dnb0=[];
ind_dl0=3; dl0=0;
[w_t,dnb_t,l0_t,Ia_new,Ib_new]=tenseg_load_prestress(substep,ind_w,w,ind_dnb,dnb0,ind_dl0,dl0,l0,b,gravity,[0;9.8;0],C,mass);

%% input file of ANSYS
ansys_input_gp_00(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),ind_w,w,ind_dnb,dnb0,'half_Tbar_ansys');
% ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),'half_Tbar_ansys');
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
% data_out=static_solver(data);        %solve equilibrium using mNewton method
data_out=static_solver2(data);        %solve equilibrium using mNewton method
% data_out{i}=equilibrium_solver_pinv(data);        %solve equilibrium using mNewton method

t_t=data_out.t_out;          %member force in every step
n_t=data_out.n_out;          %nodal coordinate in every step
N_out=data_out.N_out;
%% save output data
if savedata==1
    save (['half_Tbar_',num2str(material{1}),'_slack_',num2str(material{2})]);
end

%% plot member force
tenseg_plot_result(1:substep,t_t(1:3,:),{'element 1','element 2','element 3'},{'Load step','Force (N)'},'plot_member_force.png',saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(1:substep,n_t([3*4-2,3*4],:),{'4X','4Z'},{'Time (s)','Coordinate (m)'},'plot_coordinate.png',saveimg);

%% Plot final configuration
tenseg_plot_catenary( reshape(n_t(:,end),3,[]),C_b,C_s,[],[],[0,0],[],[],l0_t(index_s,end));

%% make video of the dynamic
name=['half_Tbar_',num2str(material{1}),'_slack_',num2str(material{2})];
% tenseg_video(n_t,C_b,C_s,[],substep,name,savevideo);
tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,[],[],[],min(substep,50),name,savevideo,material{2});

%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
