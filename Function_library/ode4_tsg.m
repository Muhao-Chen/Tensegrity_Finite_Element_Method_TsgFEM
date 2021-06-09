function data_out = ode4_tsg(odefun,tspan,y0,data_in)
%ODE4  Solve differential equations with a non-adaptive method of order 4.
%   Y = ODE4(ODEFUN,TSPAN,Y0) with TSPAN = [T1, T2, T3, ... TN] integrates
%   the system of differential equations y' = f(t,y) by stepping from T0 to
%   T1 to TN. Function ODEFUN(T,Y) must return f(t,y) in a column vector.
%   The vector Y0 is the initial conditions at T0. Each row in the solution
%   array Y corresponds to a time specified in TSPAN.
%
%   Y = ODE4(ODEFUN,TSPAN,Y0,P1,P2...) passes the additional parameters
%   P1,P2... to the derivative function as ODEFUN(T,Y,P1,P2...).
%
%   This is a non-adaptive solver. The step sequence is determined by TSPAN
%   but the derivative function ODEFUN is evaluated multiple times per step.
%   The solver implements the classical Runge-Kutta method of order 4.
%
%   Example
%         tspan = 0:0.1:20;
%         y = ode4(@vdp1,tspan,[2 0]);
%         plot(tspan,y(:,1));
%     solves the system y' = vdp1(t,y) with a constant step size of 0.1,
%     and plots the first component of the solution.
%%
global l f n f_int l_int n_d stress strain
data_out=data_in;  %initialize output data
out_tspan=data_in.out_tspan;
ne=data_in.ne;
nn= data_in.nn;
% l_int=data_in.l0_t(:,1);           % initialize the intermidiate variable
% f_int=zeros(size(data_in.l0_t(:,1)));
l_int=data_in.l0(:,1);           % initialize the intermidiate variable
f_int=zeros(size(data_in.l0(:,1)));
%%
if ~isnumeric(tspan)
    error('TSPAN should be a vector of integration steps.');
end

if ~isnumeric(y0)
    error('Y0 should be a vector of initial conditions.');
end

h = diff(tspan);
if any(sign(h(1))*h <= 0)
    error('Entries of TSPAN are not in order.')
end

try
    f0 = feval(odefun,tspan(1),y0,data_in);
catch
    msg = ['Unable to evaluate the ODEFUN at t0,y0. ',lasterr];
    error(msg);
end
y0 = y0(:);   % Make a column vector.
if ~isequal(size(y0),size(f0))
    error('Inconsistent sizes of Y0 and f(t0,y0).');
end

neq = length(y0);
N = length(tspan);
Y = zeros(neq,N);
F = zeros(neq,4);
%initialize time history data
data_out.t_t=zeros(ne,numel(out_tspan));
data_out.n_t=zeros(3*nn,numel(out_tspan));
data_out.l_t=zeros(ne,numel(out_tspan));
Y(:,1) = y0;
for i = 2:N
    ti = tspan(i-1);
    hi = h(i-1);
    yi = Y(:,i-1);
    F(:,1) = feval(odefun,ti,yi,data_in);
    if  sum(out_tspan==ti)
        disp(ti);
        data_out.t_t(:,out_tspan==ti)=f;      %member force
        data_out.n_t(:,out_tspan==ti)=n;
        data_out.l_t(:,out_tspan==ti)=l; 
        data_out.nd_t(:,out_tspan==ti)=n_d; 
        data_out.stress_t(:,out_tspan==ti)=stress; 
        data_out.strain_t(:,out_tspan==ti)=strain; 
    end
        f_int=f; l_int=l;               %store the force and length(for plastic calculation)
    F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),data_in);
    F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),data_in);
    F(:,4) = feval(odefun,tspan(i),yi+hi*F(:,3),data_in);
    Y(:,i) = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));
end
%% calculate the time history data
% %initialize time history data
% data_out.t_t=zeros(ne,numel(out_tspan));
% data_out.n_t=zeros(size(Ia,1),numel(out_tspan));
%   for i=1:numel(out_tspan)
%       %record info of structrue in step i
% %   [l,q,E,f,n,n_d]=tenseg_dyn_info(out_tspan(i),Y(:,  find(tspan==out_tspan(i))),data_in);
% [~,~,~,f,n,~]=tenseg_dyn_info(out_tspan(i),Y(:,tspan==out_tspan(i)),data_in);%record f n
% %     data_out.l_t(:,i)=l;
% %     data_out.q_t(:,i)=q;
% %     data_out.E_t(:,i)=E;
%     data_out.t_t(:,i)=f;      %member force
%     data_out.n_t(:,i)=n;
% %     data_out.n_d_t(:,i)=n_d;
%   end
if  sum(out_tspan==tspan(end))
    disp(ti);
    data_out.t_t(:,out_tspan==tspan(end))=f;      %member force
    data_out.n_t(:,out_tspan==tspan(end))=n;
    data_out.l_t(:,out_tspan==tspan(end))=l;
    data_out.stress_t(:,out_tspan==tspan(end))=stress;
    data_out.strain_t(:,out_tspan==tspan(end))=strain;
end
% feval(odefun,tspan(end),Y(:,end),data_in);
%   data_out.t_t(:,end)=f;      %member force
%     data_out.n_t(:,end)=n;
%% output data
data_out.Ya_t=Y;
% data_out.na_t=Y(1:neq/2,:);
% data_out.na_d_t=Y(neq/2+1:end,:);
end
