function gr=generate_gr_gr1(gr1,p)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function input group cell of the first unit, number of unit p,
% output gr index.
%
% Inputs:
%   gr1: group index
%	p: structure complexity
%
% Outputs:
%	gr: group matrix
%%
n_g = size(gr1,1);
m1 = 0;
for i = 1:n_g
    m1 = m1+numel(gr1{i,1});   % m1 number of members in one group
end

for i = 1:n_g
    gr{i,1} = kron(m1*[0:p-1],ones(size(gr1{i,1})))+kron(ones(1,p),gr1{i,1});
end
end