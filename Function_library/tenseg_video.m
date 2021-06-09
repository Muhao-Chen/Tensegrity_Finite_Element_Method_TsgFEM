function tenseg_video(n_t,C_b,C_s,axislim,num,name,savevideo,R3Ddata)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function make video for dynamic simulation.
%
% Inputs:
%   n_t: time history of tensegrity configuration
%   axislim: limit of axis, like [-20,20,-15,15,0,70], or use []
%   num: number of pictures in the video
%   name: name of the file
%   savevideo: save video or not
% Outputs:
%	make the video
%%
% automatic axislim
if isempty(axislim)
    %     axislim=='auto'
    nx_t=n_t(1:3:end); min_x=min(min(nx_t));max_x=max(max(nx_t));d_x=max_x-min_x+5e-2;
    ny_t=n_t(2:3:end); min_y=min(min(ny_t));max_y=max(max(ny_t));d_y=max_y-min_y+5e-2;
    nz_t=n_t(3:3:end); min_z=min(min(nz_t));max_z=max(max(nz_t));d_z=max_z-min_z+5e-1;
    
    d=max(max(n_t))-min(min(n_t));
    axislim=[min_x-0.2*d_x,max_x+0.2*d_x,min_y-0.2*d_y,max_y+0.2*d_y,min_z-0.2*d_z,max_z+0.2*d_z];
end

if nargin ==7 % no radius
    R3Ddata=[];
end

if savevideo==1
    figure(99);
    %     set(gcf,'Position',get(0,'ScreenSize')); % full screen
    for p = 1:floor(size(n_t,2)/num):size(n_t,2)
        N = reshape(n_t(:,p),3,[]);
        %         tenseg_plot(N,C_b,C_s,99,[],[],[],R3Ddata);hold on
        %         delete_index = [7,13];
        %         [N,C_b,C_s] = tenseg_delete_extra_nodes(delete_index,N,C_b,C_s);
        tenseg_plot(N,C_b,C_s,99,[],[],[],R3Ddata);hold on
        set(gcf,'color','w');
        axis(axislim)
        tenseg_savegif_forever(name);
        hold off;
    end
    close
end
end
