%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%An double layer tensegrity tower with simplex%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Wood','Rubber_band');
material{1}='linear_elastic'; % index for material properties: multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
dt=0.001;               % time step in dynamic simulation
auto_dt=1;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=1;                   % final time of dynamic simulation
out_dt=0.001;            % output data interval(approximately, not exatly)
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Specify node positions
N = [-1 0 0; 0 -1 0; 1 0 0; 0 1 0]';

% Specify bar connectivity
Cb_in = [1 3; 2 4];
C_b = tenseg_ind2C(Cb_in,N);

% Specify string connectivity
Cs_in = [1 2; 2 3; 3 4;4 1];  % String one is node 1 to 2
C_s = tenseg_ind2C(Cs_in,N);
C=[C_b;C_s];
[ne,nn]=size(C);        % ne:No.of element;nn:No.of node

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
title('Tbar');

%% Boundary constraints
pinned_X=([])'; pinned_Y=([])'; pinned_Z=([])';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% Group information
%generate group index
gr={[1 2];[3 5];[4 6]};     % number of elements in one group
% gr=[];                     %if no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% self-stress design
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%SVD of equilibrium matrix
[U1,U2,V1,V2,S1]=tenseg_svd(A_1ag);

%external force in equilibrium design
w0=zeros(numel(N),1); w0a=Ia'*w0;

%prestress design
% index_gp=[1,2];                 % number of groups with designed force
% fd=[-1e5,1e4];              % force in bar is given as -1000
% [q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design
q_gp=[-1e3;1e3;1e3];
t_gp=l_gp.*q_gp;
q=Gp*q_gp;
t=Gp*t_gp;

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
% A_gp=[1e-3;1e-3;1e-3];
% A=Gp*A_gp;
% r_b=0.5;
% r_s=0.1;

[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
l0([3,5],1)=0.7*l0([3,5],1);
% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.01,0.03],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.01,0.03],r_s);
R3Ddata.Nradius=0.04*ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Tbar',R3Ddata);

%% input file of ANSYS
ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),fullfile(savePath,'Tbar'));

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% mode analysis
num_plt=4:6;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);

%% input data for dynamic analysis
% time step
if auto_dt
dt=pi/(64*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
% dt=dt*1e-1;
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

[w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,b,'step',gravity,[0;0;0],C,mass,3*3-1,0);

% give initial speed of free coordinates
n0a_d=2*ones(numel(a),1);        %initial speed in X direction

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
% tenseg_plot_exforce(Ib,tspan,w_t,[1],dnb_t,dnb_d_t,dnb_dd_t,[1],saveimg);

%% dynamic analysis 
% solve dynamic equation
data_out=dynamic_solver(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate 
l_t=data_out.l_t;   %time history of members' length 
nd_t=data_out.nd_t;   %time history of nodal coordinate

%% save output data
if savedata==1
    save (['Tbar_dynamics',material{1},'.mat']);
end

%% plot member force 
tenseg_plot_result(out_tspan,t_t(1:6,:),{'element 1','element 2','element 3','element 4','element 5','element 6'},...
                   {'Time (s)','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t(3*1-2,:),{'1X'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*1-1,:),{'1Y'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Y_coordinate.png'),saveimg);

%% make video of the dynamic
name=fullfile(savePath,['Tbar','tf_',num2str(tf),material{1}]);
tenseg_video(n_t,C_b,C_s,[],100,name,savevideo);

%% linearized dynaimcs 
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);
