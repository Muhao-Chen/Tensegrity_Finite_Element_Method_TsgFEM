function tenseg_savegif_forever(varargin)
%% tenseg_savegif_forever
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function will save the animation in the for loop as gif file.
% How to use it?
% 1. When use it, it has to be in for loop, the default name is untitledgif.gif?;
% 2. tenseg_savegif_forever('filename') will named the animation as filename.gif?.
% 3. tenseg_savegif_forever('filename.gif') will named the animation as filename.gif?.
% Note: please use 'clear all' to clean all the global variables when using this function.
%%
% warning off;

%Parse input 
switch nargin
    case 1
        filename=varargin{1};
        dt=0.1; %default dt=0.1
    case 2
        filename=varargin{1};
        dt=varargin{2};
end

loops=65535;
global iiii; global time;
p=clock;

if isempty(filename)
    filename='untitledgif.gif';
else    
    if length(filename)<4
        filename=[filename,'.gif'];
    else
        if ~strcmp(filename(end-3:end),'.gif')
            filename=[filename,'.gif'];
        end
    end
end
if isempty(iiii)
    iiii=0; time=p(6);
else
    iiii=iiii+1;
    if ((p(6)>=time)*(p(6)-time)+(p(6)<time)*(p(6)+60-time))>10
        iiii=0;
    end
    time=p(6);
end
f=getframe(gcf); % one may change this into gca
f=frame2im(f);
[f,map]=rgb2ind(f,256);
if iiii==0
    imwrite(f,map,filename,'LoopCount',loops,'delaytime',dt);
else
    imwrite(f,map,filename,'writemode','append','delaytime',dt);
end