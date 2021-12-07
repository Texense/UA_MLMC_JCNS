% MC_Vseries_for_network_ninput_I_based_norms.m

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

%% main function to evolve recurrent network
function [Vs,spike,t,rates,V_Nm] = MC_Vseries_for_network_ninput_I_based_norms(dt,T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  g_L = 20;

  

  N = length(V_0s);
% Read input spike time points in eventsN, and put the points in numerical
% iteration process
  [A,B,Times] = find(eventsN);
  R = sparse(A,floor(Times/dt)+1,ones(size(A)),N,length(t));

  % add refractory period
  t_ref = 1;
  N_ref = floor(t_ref/dt); % number of steps in ref period
  spike = zeros(N,length(t));
 
  % matrix Vs records V for every neuron on every time step
  Vs = zeros(N,length(t)+N_ref);
  Vs(:,1) = V_0s;
  Ref = ones(size(Vs)); 
  NN = length(V_0s);
  SE = [see*S(1:NE,1:NE),zeros(size(S(1:NE,NE+1:NN)));sie*S(NE+1:NN,1:NE),zeros(size(S(NE+1:NN,NE+1:NN)))];
  SI = [zeros(size(S(1:NE,1:NE))),sei*S(1:NE,NE+1:NN);zeros(size(S(NE+1:NN,1:NE))),sii*S(NE+1:NN,NE+1:NN)]; 
  
  for ii = 1:steps
   
%    dV = (dt*(-Vs(:,ii)/g_L) +(fs.*R(:,ii)+SE*spike(:,ii)-SI*spike(:,ii))/g_L ).*Ref(:,ii);
   dV = (dt*(-Vs(:,ii)/g_L) +(fs.*R(:,ii)+SE*spike(:,ii)-SI*spike(:,ii))).*Ref(:,ii);
   VV = Vs(:,ii) + dV; % the rough output
   V_pls = mod(subplus(VV),1);
   spike(:,ii+1) = subplus(VV)-V_pls;
   Ref(:,ii+1:ii+N_ref) = Ref(:,ii+1:ii+N_ref) - repmat(spike(:,ii+1),[1,N_ref]);
   Vs(:,ii+1) = min(V_pls,VV);
   
  end
  % use norm_rate function to calculate firing rate m(t) for the whole network 
  pieces = floor(T/500); %500ms is the time length for each piece.
  %As in norm_rate function each entry in "rates" cover consecutive 2 pieces, actually the unit
  %for rates is # of spikes for 600 neurons per second.
  rates = SpikeCount_slide(t,spike,5,1);%norm_rate(t,spike,pieces);
  % WinSize 5ms, slide 1ms
  V_Nm = norm_voltage(t,Vs,V_0s,pieces);  
  
end