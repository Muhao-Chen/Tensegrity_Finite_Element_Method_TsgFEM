function tenseg_video_slack(n_t,C_b,C_s,l0_t,index_s,R3Ddata,view_vec,axislim,num,name,savevideo,slack)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function make video for dynamic simulation(plot slack string as a
% catenary.
%
% Inputs:
%   n_t: time history of tensegrity configuration
%   C_b,C_s: connectivity matrix of bar and string
%   axislim: limit of axis, like [-20,20,-15,15,0,70], or use [] 
%   num: number of pictures in the video
%   name: name of the file
%   savevideo: save video or not
%   slack: plot strings slack as cantenary
% Outputs:
%	make the video
%%
% automatic axislim
if isempty(axislim)
%     axislim=='auto'
nx_t=n_t(1:3:end); min_x=min(min(nx_t));max_x=max(max(nx_t));d_x=max_x-min_x+1e-3;
ny_t=n_t(2:3:end); min_y=min(min(ny_t));max_y=max(max(ny_t));d_y=max_y-min_y+1e-3;
nz_t=n_t(3:3:end); min_z=min(min(nz_t));max_z=max(max(nz_t));d_z=max_z-min_z+1e-3;

d=max(max(n_t))-min(min(n_t));
axislim=[min_x-0.02*d_x,max_x+0.02*d_x,min_y-0.02*d_y,max_y+0.02*d_y,min_z-0.02*d_z,max_z+0.02*d_z];
end

if savevideo==1
    figure(99);
    set(gcf,'Position',get(0,'ScreenSize'));          %full screen
    for p = 1:floor(size(n_t,2)/num):size(n_t,2)
        switch slack
            case 0
        tenseg_plot(reshape(n_t(:,p),3,[]),C_b,C_s,99,[],view_vec,[],R3Ddata);hold on
            case 1
 tenseg_plot_catenary( reshape(n_t(:,p),3,[]),C_b,C_s,99,[],view_vec,[],R3Ddata,l0_t(index_s,p));hold on
        end
        set(gcf,'color','w');
        axis(axislim)
        tenseg_savegif_forever(name);
        hold off;
    end
    close 
end
end

