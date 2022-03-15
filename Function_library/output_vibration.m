function output_vibration(data,name)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the time history of vibration of ground motion to a
% *.txt file
%
% Inputs:
%   data: a column vector containing acceleration of ground motion
%   name: a string containing the name of output txt file
% Outputs:
%   a.txt file
%%
fid11=fopen(name,'w');
for i=1:numel(data)
    fprintf(fid11,'%25.25f \n',data(i));
end
fclose(fid11);
end

