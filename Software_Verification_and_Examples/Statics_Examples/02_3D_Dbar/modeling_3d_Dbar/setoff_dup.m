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
index = [];
for i=1:length(N(:,1))
    a = find(((N(:,1)-N(i,1)).^2 + (N(:,2)-N(i,2)).^2 +(N(:,3)-N(i,3)).^2).^0.5<10^-precision ==1)
    tmpa = a(find(a>i));
    if ~isempty(tmpa)
    index = [index;tmpa];
    end
end
N(index,:)=[];
end