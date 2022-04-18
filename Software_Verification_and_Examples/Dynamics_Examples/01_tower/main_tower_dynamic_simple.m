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
%%
%EXAMPLE
clc;clear;close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Steel_Q345','Steel_string');
material{1}='linear_elastic'; % index for material properties: linear_elastic, multielastic, plastic.
material{2}=1; % index for considering slack of string (1) for yes,(0) for no (for comparision with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
amplitude=50;            % amplitude of external force of ground motion
period=0.5;             % period of seismic signal

dt=0.001;               % time step in dynamic simulation
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=0.02;            % output data interval(approximately, not exatly)
lumped=0;               % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=1;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0)
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of double layer prism.
R=10; h=30; p=3;        % radius; height; number of edge
beta=180*(0.5-1/p); 	% rotation angle

% for i=1:p               % nodal coordinate matrix N
%     N(:,i)=R*[cos(2*pi*(i-1)/p),sin(2*pi*(i-1)/p),0];
% end
% for i=p+1:2*p
%     N(:,i)=[R*cos(2*pi*(i-1)/p+beta*pi/180),R*sin(2*pi*(i-1)/p+beta*pi/180),h];
% end
% for i=2*p+1:3*p
%     N(:,i)=[R*cos(2*pi*(i-1)/p+2*beta*pi/180),R*sin(2*pi*(i-1)/p+2*beta*pi/180),2*h];
% end
angle1=2*pi*((1:p)-1)./p;
N=R*[cos(angle1); sin(angle1); zeros(1,p)];
angle2=2*pi*((1:p)-1)./p+beta*pi/180;
N=[N,[R*[cos(angle2); sin(angle2)]; h*ones(1,p)]];
angle3=2*pi*((1:p)-1)./p+2*beta*pi/180;
N=[N,[R*[cos(angle3); sin(angle3)]; 2*h*ones(1,p)]];
% Manually specify connectivity indices.
C_s_in = [4 5;5 6;6 4;7 8;8 9;9 7;1 4;2 5;3 6;4 7;5 8;6 9];  % This is indicating that string connection
C_b_in = [1 5;2 6;3 4;5 9;6 7;4 8];  % Similarly, this is saying bar 1 connects node 1 to node 2,

% Convert the above matrices into full connectivity matrices.
C_b = tenseg_ind2C(C_b_in,N);%%
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
index_gp=[1,2];                 % number of groups with designed force
fd=-1e5*ones(2,1);              % force in bar is given as -1000
[q_gp,t_gp,q,t]=tenseg_prestress_design(Gp,l,l_gp,A_1ag,V2,w0a,index_gp,fd);    %prestress design

%% cross sectional design
index_b=find(t<0);              % index of bar in compression
index_s=setdiff(1:ne,index_b);	% index of strings
[A_b,A_s,A_gp,A,r_b,r_s,r_gp,radius,E,l0,rho,mass]=tenseg_minimass(t,l,Gp,sigmas,sigmab,Eb,Es,index_b,index_s,c_b,c_s,rho_b,rho_s,thick,hollow_solid);
% Plot the structure with radius
R3Ddata.Bradius=interp1([min(radius),max(radius)],[0.3,1],r_b);
R3Ddata.Sradius=interp1([min(radius),max(radius)],[0.3,1],r_s);
R3Ddata.Nradius=ones(nn,1);
tenseg_plot(N,C_b,C_s,[],[],[],'Double layer prism',R3Ddata);

%% input file of ANSYS
ansys_input_gp(N,C,A_gp,t_gp,b,Eb,Es,rho_b,rho_s,Gp,index_s,find(t_gp>0),fullfile(savePath,'tower'));

%% mass matrix and damping matrix
M=tenseg_mass_matrix(mass,C,lumped); % generate mass matrix
% damping matrix
d=0;     %damping coefficient
D=d*2*max(sqrt(mass.*E.*A./l0))*eye(3*nn);    %critical damping

%% free vibration mode analysis
num_plt=1:3;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg,100);

%% input data for dynamic analysis
% time step
if auto_dt
    dt=pi/(8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span
%% External force and boundary constraints
% calculate external force and forced motion of nodes
[w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,b,'vib_force',gravity,[0;0;9.8],C,mass,3*(4:9)-2,50,period);
% give initial speed of free coordinates
n0a_d=0*kron(ones(numel(a)/3,1),[1;0;0]);        %initial speed in X direction

% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'impluse',gravity,[0;0;9.8],C,mass,3*(4:9)-1,1e5,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'step',gravity,[0;0;9.8],C,mass,3*(4:9)-2,1e5,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'ramp',gravity,[0;0;9.8],C,mass,3*(4:9),-1e6,period);
% [w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,a,b,'vib_nodes',gravity,[0;0;9.8],C,mass,3*(1:3)-2,50,period);

%% input data
data.N=N; data.C=C; data.ne=ne; data.nn=nn; data.Ia=Ia; data.Ib=Ib;
data.E=E; data.A=A; data.l0=l0; data.index_b=index_b; data.index_s=index_s;
data.consti_data=consti_data;   data.material=material; %constitue info
data.w_t=w_t;  % external force
data.dnb_t=dnb_t; data.dnb_d_t=dnb_d_t;  data.dnb_dd_t=dnb_dd_t; % forced movement of pinned nodes
data.n0a_d=n0a_d;        %initial speed of free coordinates
data.M=M;data.D=D;
data.tf=tf;data.dt=dt;data.tspan=tspan;data.out_tspan=out_tspan;

%% output the seismic data
output_vibration(dz_a_t(1,:)',fullfile(savePath,'tjx.txt'));

%% plot external force information
% tenseg_plot_exforce(Ib,tspan,w_t,(4:6),dnb_t,dnb_d_t,dnb_dd_t,(1),saveimg);

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate
l_t=data_out.l_t;   %time history of members' length
nd_t=data_out.nd_t;   %time history of nodal velocity

%% save output data
if savedata==1
    save (fullfile(savePath,['FEM_tower_simple',material{1},'.mat']));
end

%% plot member force
tenseg_plot_result(out_tspan,t_t(7:8,:),{'element 7','element 8'},{'Time (s)','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);

%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t(3*8-2,:),{'8X'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*8-1,:),{'8Y'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Y_coordinate.png'),saveimg);
%% plot velocity of nodal coordinate
tenseg_plot_result(out_tspan(1:end-1),nd_t(3*8-1,:),{'8Y'},{'Time (s)','Velocity (m/s)'},fullfile(savePath,'Y_coordinate.png'),saveimg);
%% plot member length
tenseg_plot_result(out_tspan,l_t(7:8,:),{'element 7','element 8'},{'Time (s)','Length (m)'},fullfile(savePath,'member_length.png'),saveimg);
%% make video of the dynamic
name=fullfile(savePath,['tower','tf_',num2str(tf),material{1},num2str(material{2})]);
% tenseg_video(n_t,C_b,C_s,[],50,name,savevideo,R3Ddata);
tenseg_video(n_t,C_b,C_s,[],50,name,savevideo);
%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);