function N = setoff_dup(N,precision)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% This function is going to delete duplicate nodes
%  Inputs:
%       N: N*3 matrix
%       precision: Integer for precision
%  Outputs:
%       N: N*3 matrix with no duplicates
%%

[~,indUni]=unique(round(N,precision),'rows'); %Indices of unique points
N=N(indUni,:); %Keep only unique points

end