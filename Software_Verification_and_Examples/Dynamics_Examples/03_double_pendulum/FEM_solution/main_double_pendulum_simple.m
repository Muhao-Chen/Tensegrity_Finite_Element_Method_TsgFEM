%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Double pendulum%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% [1] structure design(calculate equilibrium matrix,
% group matrix,prestress mode, minimal mass design)
% [2] modal analysis(calculate tangent stiffness matrix, material
% stiffness, geometry stiffness, generalized eigenvalue analysis)
% [3] dynamic simulation

%EXAMPLE
clc;clear;close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

auto_dt=1;              % use(1 or 0) auto time step, converengency is guaranteed if used
out_dt=1e-4;            % output data interval(approximately, not exatly)
tf=5;                   % final time of dynamic simulation
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=1;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of structure
% initial configuration of structure
n1=[0,0,0]';
n2=[sqrt(2)/2,-sqrt(2)/2,0]';
n3=[sqrt(2),-sqrt(2),0]';
N=[n1,n2,n3];

%connectivity of structure
C_b_in=[1,2;2,3];
C_s_in=[];
% Convert above index notation into actual connectivity matrices
C_b = tenseg_ind2C(C_b_in,N);
C_s = tenseg_ind2C(C_s_in,N);
C=[C_b;C_s];

[ne,nn]=size(C);% ne:No.of element;nn:No.of node
% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);

%% Boundary constraints
pinned_X=(1)'; pinned_Y=(1)'; pinned_Z=(1:3)';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];                      %no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% equilibrium matrix
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%% Cross sectional & young's modulus
index_b=1:ne; index_s=[];       % define bar, string index
A=1e-4*ones(ne,1);            % define cross sectional area vector
E=Eb*ones(ne,1);                % define Young's modulus vector

%%  members' force & rest length
t=0*ones(ne,1);                 % zero internal force
l0=E.*A.*l./(t+E.*A);

%% input data for dynamic analysis
% density vector
rho=rho_b*ones(ne,1);
% mass matrix
mass=rho.*A.*l0;           % mass vector
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=1:3;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% input data for dynamic analysis
% time step
if auto_dt
    dt=pi/(16*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_dt=dt;            % output data interval(approximately, not exatly)
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and forced motion of nodes
[w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,b,'step',gravity,[0;9.8;0],C,mass,3*3-1,0);

% give initial speed of free coordinates
n0a_d=0*ones(numel(a),1);        %initial speed in X direction

% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% plot external force information
% tenseg_plot_exforce(Ib,tspan,w_t,[2,3],dnb_t,dnb_d_t,dnb_dd_t,[1],saveimg);

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate
l_t=data_out.l_t;   %time history of members' length
%% save output data
if savedata==1
    save (fullfile(savePath,['FEM_double_pendulum','.mat']));
end

%% plot member force
tenseg_plot_result(out_tspan,t_t(1:2,:),{'element 1','element 2'},{'Time (s)','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t(3*[2:3]-2,:),{'2X','3X'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*[2:3]-1,:),{'2Y','3Y'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Y_coordinate.png'),saveimg);

%% Plot error of bar length
tenseg_plot_result(out_tspan,l_t(1:2,:)-1,{'Bar 1','Bar 2'},{'Time (s)','Error of Bar Length (m)'},fullfile(savePath,'Error of bar length.png'),saveimg);

%% make video of the dynamic
name=fullfile(savePath,['double_pendulum','tf_',num2str(tf),material{1}]);
tenseg_video(n_t,C_b,C_s,[-1.5,1.5,-2,0.1,-1,1],100,name,savevideo);

%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
