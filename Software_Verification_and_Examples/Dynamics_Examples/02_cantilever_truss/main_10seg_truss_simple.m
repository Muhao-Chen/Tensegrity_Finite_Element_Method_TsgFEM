% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Cantilever truss%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% [1] structure design(give N, C, E, A, rho, l0)
% [2] modal analysis(calculate tangent stiffness matrix, material
% stiffness, geometry stiffness, generalized eigenvalue analysis
% to get frequency and mode shape)
% [3] dynamic simulation in step load

%EXAMPLE
clc; clear; close all;
% Global variable
[consti_data,Eb,Es,sigmab,sigmas,rho_b,rho_s]=material_lib('Q345_blin','Steel_string');
material{1}='plastic'; % index for material properties: linear_elastic, multielastic, plastic.
material{2}=0; % index for considering slack of string (1) for yes,(0) for no (for compare with ANSYS)

% cross section design cofficient
thick=6e-3;        % thickness of hollow bar
hollow_solid=0;          % use hollow bar or solid bar in minimal mass design (1)hollow (0)solid
c_b=0.1;           % coefficient of safty of bars 0.5
c_s=0.1;           % coefficient of safty of strings 0.3

% dynamic analysis set
amplitude=-1e5;            % amplitude of external force of ground motion

dt=1e-4;                % artifitial given time step
auto_dt=0;              % use(1 or 0) auto time step, converengency is guaranteed if used
tf=2;                   % final time of dynamic simulation
out_dt=1e-2;            % output data interval(approximately, not exatly)
lumped=0;    % use lumped matrix 1-yes,0-no
saveimg=0;              % save image or not (1) yes (0)no
savedata=0;             % save data or not (1) yes (0)no
savevideo=0;            % make video(1) or not(0)
gravity=0;              % consider gravity 1 for yes, 0 for no
% move_ground=0;          % for earthquake, use pinned nodes motion(1) or add inertia force in free node(0)
savePath=fullfile(fileparts(mfilename('fullpath')),'data_temp'); %Save files in same folder as this code

%% N C of the structure
% Manually specify node positions of double layer prism.
N=[0 0 0;1 0 0;2 0 0;3 0 0;4 0 0;5 0 0;6 0 0;7 0 0;8 0 0;9 0 0;10 0 0;
    0 1 0;1 1 0;2 1 0;3 1 0;4 1 0;5 1 0;6 1 0;7 1 0;8 1 0;9 1 0;10 1 0]';   %nodal coordinate
C_in=[1 2; 2 3;3 4;4 5;5 6;6 7;7 8;8 9;9 10;10 11;
    1 12;1 13;2 13;2 14;3 14;3 15;4 15;4 16;5 16;5 17;6 17;6 18;7 18;7 19;8 19;8 20;9 20;9 21;10 21;10 22;11 22;
    12 13;13 14;14 15;15 16;16 17;17 18;18 19;19 20;20 21;21 22];
C = tenseg_ind2C(C_in,N);

[ne,nn]=size(C);        % ne:No.of element;nn:No.of node
C_b=C;C_s=[];

% Plot the structure to make sure it looks right
tenseg_plot(N,C_b,C_s);
% title('cantilever truss');
%% Boundary constraints
pinned_X=[1,12]'; pinned_Y=[1,12]'; pinned_Z=[1:22]';
[Ia,Ib,a,b]=tenseg_boundary(pinned_X,pinned_Y,pinned_Z,nn);

%% generate group index for tensegrity torus structure
gr=[];                      %no group is used
Gp=tenseg_str_gp(gr,C);    %generate group matrix

%% equilibrium matrix
%Calculate equilibrium matrix and member length
[A_1a,A_1ag,A_2a,A_2ag,l,l_gp]=tenseg_equilibrium_matrix1(N,C,Gp,Ia);

%% Cross sectional & young's modulus
index_b=1:ne; index_s=[];       % define bar, string index
A=0.0025*ones(ne,1);            % define cross sectional area vector
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

%% output for anlysis
ansys_input_gp(N,C,A,t,b,Eb,Es,rho_b,rho_s,Gp,index_s,index_s,fullfile(savePath,'truss_10seg2.txt'));

%% mode analysis
num_plt=1:3;        % number of modes to plot
[V_mode,omega]=tenseg_mode(Ia,C,C_b,C_s,N(:),E,A,l0,M,num_plt,saveimg);
%% input data for dynamic analysis
% time step
if auto_dt
    dt=pi/(8*max(omega)); 	% time step dt is 1/8 of the smallest period, guarantee convergence in solving ODE
end
tspan=0:dt:tf;
out_tspan=interp1(tspan,tspan,0:out_dt:tf, 'nearest','extrap');  % output data time span

% calculate external force and forced motion of nodes
[w_t,dnb_t,dnb_d_t,dnb_dd_t,dz_a_t]=tenseg_ex_force(tspan,b,'step',gravity,[0;0;9.8],C,mass,3*22-1,amplitude);

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
% tenseg_plot_exforce(Ib,tspan,w_t,22,dnb_t,dnb_d_t,dnb_dd_t,1,saveimg);

%% dynamic analysis
% solve dynamic equation
data_out=dynamic_solver(data);        %solve ODE of dynamic equation
% time history of structure
t_t=data_out.t_t;   %time history of members' force
n_t=data_out.n_t;   %time history of nodal coordinate
l_t=data_out.l_t;   %time history of members' length
stress_t=data_out.stress_t; % time history of members' stress
strain_t=data_out.strain_t; % time history of members' strain

%% save output data
if savedata==1
    save (fullfile(savePath,['FEM_cantilever_10seg',material{1},'.mat'])) ;
end

%% plot member force
tenseg_plot_result(out_tspan,t_t(1:2,:),{'element 1','element 2'},{'Time (s)','Force (N)'},fullfile(savePath,'member_force.png'),saveimg);
%% video member stress-strain
tenseg_video_strain_stress(data_out,[1,32,37],0.1);

%% plot member stress-strain
% time=[1/3,2/3,1];
% time=0.5
time=1;
tenseg_plot_strain_stress(data_out,[1,32,37],time);
%plot stress
tenseg_plot_result(out_tspan,stress_t([1,32,37],:),{'1','32','37'},{'Time (s)','Stress (Pa)'},fullfile(savePath,'Stress.png'),saveimg);


%% Plot nodal coordinate curve X Y
tenseg_plot_result(out_tspan,n_t(3*22-2,:),{'22X'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'X_coordinate.png'),saveimg);
tenseg_plot_result(out_tspan,n_t(3*22-1,:),{'22Y'},{'Time (s)','Coordinate (m)'},fullfile(savePath,'Y_coordinate.png'),saveimg);

%% Plot configuration
%plot configuration in one given time
[~,iii]=min(abs(out_tspan-1));       %this is to selcet the No. related to 0.35s
tenseg_plot(reshape(n_t(:,iii),3,[]),C_b,C_s,[],[],[0,90]);   %plot structure in 0.3s
axis([0,10.5,-1.5,1.2,0-1,1]);
%plot configuration in three given time
time=[1/3,2/3,1];
picture=figure;
pic=zeros(1,numel(time));
for i=1:numel(time)
    [~,pic(i)]=min(abs(out_tspan-time(i)));       %this is to selcet the No. related to 0s
    tenseg_plot(reshape(n_t(:,pic(i)),3,[]),C_b,C_s,picture,[],[0,90]);   %plot structure in 0.3s
    hold on;
end
axis([0,10.5,-1.5,1.2,0-1,1]);

%% make video of the dynamic
name=fullfile(savePath,['cantilever_','tf_',num2str(tf),material{1}]);
tenseg_video(n_t,C_b,C_s,[0,10.5,-1.5,1.2,0-1,1],100,name,savevideo);

%% linearized dynaimcs
[A_lin,B_lin]=tenseg_lin_mtrx(C,N(:),E,A,l0,M,D,Ia,A_1a);