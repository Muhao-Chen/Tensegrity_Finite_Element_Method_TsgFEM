function d=eulerdst(X,Y)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
sumd=0;
for i=1:length(X)
    sumd = sumd + (X(i)-Y(i))^2;
end
d = sqrt(sumd);
end