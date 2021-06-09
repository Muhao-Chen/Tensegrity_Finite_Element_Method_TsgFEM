function [x,y,z] = catenary_3D(A,B,r_length,N,sagInit)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% given two points A=[ax ay az] and B=[bx by bz] in space,
% rope length r_length, and the number of intermediate points N,
% outputs the coordinates x ,y and z of the hanging rope from a to b
% the optional input sagInit initializes the sag parameter for the
% root-finding procedure.
%%
maxIter	= 100;       % maximum number of iterations
minGrad	= 1e-10;     % minimum norm of gradient
minVal	= 1e-8;      % minimum norm of sag function
stepDec	= 0.5;       % factor for decreasing stepsize
minStep	= 1e-9;		 % minimum step size
minHoriz	= 1e-3;		 % minumum horizontal distance

%% transform a 3D problem to 2D problem
a=[0,A(3)];
b=[norm(A(1:2)-B(1:2)),B(3)];

%%
if nargin < 5
	sag = 1;
else
	sag = sagInit;
end

if a(1) > b(1)
	[a,b]	= deal(b,a);
end

d = b(1)-a(1);
h = b(2)-a(2);


if abs(d) < minHoriz % almost perfectly vertical
	X = ones(1,N)*(a(1)+b(1))/2;
	if r_length < abs(h) % rope is stretched
		Y		= linspace(a(2),b(2),N);
	else
		sag	= (r_length-abs(h))/2;
		n_sag = ceil( N * sag/r_length );
		y_max = max(a(2),b(2));
		y_min = min(a(2),b(2));
		Y		= linspace(y_max,y_min-sag,N-n_sag);
		Y		= [Y linspace(y_min-sag,y_min,n_sag)];
    end
    %% transform from 2D to 3D   
    x=ones(1,N)*A(1);
    y=ones(1,N)*A(2);
    z=Y;
	return;
end

X = linspace(a(1),b(1),N);

if r_length <= sqrt(d^2+h^2) % rope is stretched: straight line
	Y = linspace(a(2),b(2),N);
	    %% transform from 2D to 3D   
    x=ones(1,N)*A(1)+(B(1)-A(1))/norm(A(1:2)-B(1:2))*X;
    y=ones(1,N)*A(2)+(B(2)-A(2))/norm(A(1:2)-B(1:2))*X;
    z=Y;
else
	% find rope sag
	g  = @(s) 2*sinh(s*d/2)/s - sqrt(r_length^2-h^2);
	dg = @(s) 2*cosh(s*d/2)*d/(2*s) - 2*sinh(s*d/2)/(s^2);

	for iter = 1:maxIter
		val		= g(sag); 
		grad		= dg(sag);
		if abs(val) < minVal || abs(grad) < minGrad
			break
		end
		search	= -g(sag)/dg(sag);
		
		alpha		= 1;
		sag_new  = sag + alpha*search;
		
		while sag_new < 0 || abs(g(sag_new)) > abs(val)
			alpha		= stepDec*alpha;
			if alpha < minStep
				break;
			end
			sag_new	= sag + alpha*search;			
		end
		
		sag = sag_new;
	end

	% get location of rope minimum and vertical bias
	x_left	= 1/2*(log((r_length+h)/(r_length-h))/sag-d);
	x_min		= a(1) - x_left;
	bias		= a(2) - cosh(x_left*sag)/sag;

	Y			= cosh((X-x_min)*sag)/sag + bias;
    
    %% transform from 2D to 3D   
    x=ones(1,N)*A(1)+(B(1)-A(1))/norm(A(1:2)-B(1:2))*X;
    y=ones(1,N)*A(2)+(B(2)-A(2))/norm(A(1:2)-B(1:2))*X;
    z=Y;
end
end