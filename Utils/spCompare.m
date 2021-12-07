% spCompare.m

% This file is part of UA_MLMC_JCNS, a collection of
% numerical experiments on multilevel Monte Carlo for
% spiking neuron networks.   It accompanies the
% paper "Multilevel Monte Carlo for Cortical Circuit Models,"
% to appear in the Journal of Computational Neuroscience.
% 
% Copyright (C) 2021 by Zhuo-Cheng Xiao
% <zx555@nyu.edu>
% 
% This program is free software; you can redistribute
% it and/or modify it under the terms of the GNU
% General Public License as published by the Free
% Software Foundation; either version 2 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it
% will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General
% Public License along with this program; if not, write
% to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA 02110-1301 USA.

%% spike comparison of two 1D spike-time series
% Assume spike2 is the more accurate one, treat it as template.
% spike1,2 should be both row vectors
% ZCX 03/11/21
function SpDiff = spCompare(spike12, spike22)
l1 = length(spike12);
l2 = length(spike22);
if l1 == l2 
    SpDiff = spike12 - spike22;
elseif max(l1,l2)>10000 % too large, just use simple solutions
    Sp1Short = l1(1:min(l1,l2));
    Sp2Short = l2(1:min(l1,l2));
    SpDiff = Sp1Short - Sp2Short;
else
    SpDiffMat = repmat(spike12',1,l2) - repmat(spike22,l1,1);
    [~,CompInd] = min(abs(SpDiffMat));
    SpDiff = zeros(size(spike22));
    for SpInd = 1:l2
        SpDiff(SpInd) = SpDiffMat(CompInd(SpInd),SpInd);
    end
end

end
