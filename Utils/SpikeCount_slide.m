% SpikeCount_slide.m

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

function rates = SpikeCount_slide(t,spike,WinSize,slideSize)
  pieces = floor((t(end)-WinSize)/slideSize);
  NWinBin = floor(WinSize/(t(2)-t(1)));
  SlideBin = floor(slideSize/(t(2)-t(1)));
  rates = zeros(size(spike,1),pieces);
  for WinInd = 1:pieces
      rates(:,WinInd) = sum(spike(:, 1+(WinInd-1)*SlideBin : (WinInd-1)*SlideBin + NWinBin),2);
  end
end