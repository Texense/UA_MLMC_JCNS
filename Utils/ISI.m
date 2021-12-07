% ISI.m

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

function [DSE,DSI] = ISI(spike, NE,NI, dt)
DSE = []; DSI = [];
for ii = 1:NE
    Sp_tii = find(spike(ii,:) == 1)*dt;
    dif_ii = Sp_tii(2:length(Sp_tii)) - Sp_tii(1:length(Sp_tii)-1);
    DSE = [DSE,dif_ii];
    clear Sp_tii dif_ii
end

for jj = (NE+1):(NE+NI)
    Sp_tjj = find(spike(jj,:) == 1)*dt;
    dif_jj = Sp_tjj(2:length(Sp_tjj)) - Sp_tjj(1:length(Sp_tjj)-1);
    DSI = [DSI,dif_jj];
    clear Sp_tjj dif_jj
end

end