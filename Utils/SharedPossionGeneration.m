% SharedPossionGeneration.m

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

%% Generate N Possion Spike trains with shared input

function eventsN = SharedPossionGeneration(lambda1,ShareFactor,T,N)

%lambda1 = 0.55;
%ShareFactor = 0.0;
Nevent    = floor(T*lambda1*2);
NSHAEvent = floor(ShareFactor*Nevent); % number of shared input 
%NIIDEvent = Nevent - NSHAEvent;        % number of iid events

eventsN = zeros(N,Nevent);
SharedEvent = Poisson_Process_2(lambda1*ShareFactor,T);
for ii = 1:N
IIDevent = Poisson_Process_2(lambda1*(1-ShareFactor),T);
FinalEvent = sort([SharedEvent;IIDevent]);
eventsN(ii,1:length(FinalEvent)) = FinalEvent';
end
eventsN(:,1) = []; % get rid of the first column, all zero.
end