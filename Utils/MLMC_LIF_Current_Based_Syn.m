% MLMC_LIF_Current_Based_Syn.m

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
function [Vs,I_Es,I_Is,spike,t,SpikeCount,V_Nm,SynIndex] = MLMC_LIF_Current_Based_Syn(dt,T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI,...
                                                                       I_E0,I_I0)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  g_L = 20;
  
  tau_E = 3;
  tau_I = 4;
  

  N = length(V_0s);
% Read input spike time points in eventsN, and put the points in numerical
% iteration process
  [A,~,Times] = find(eventsN);
  R = sparse(A,floor(Times/dt)+1,ones(size(A)),N,length(t));
  %R = full(R);
  % add refractory period
  t_ref = 1;
 
  % matrix Vs records V for every neuron on every time step
  Vs = zeros(N,length(t));     Vs(:,1)   = V_0s;
  I_Es = zeros(N,length(t));   I_Es(:,1) = I_E0;
  I_Is = zeros(N,length(t));   I_Is(:,1) = I_I0;
  %Ref = ones(size(Vs)); 
  NN = length(V_0s);
  SE = [see*S(1:NE,1:NE),zeros(size(S(1:NE,NE+1:NN)));sie*S(NE+1:NN,1:NE),zeros(size(S(NE+1:NN,NE+1:NN)))];
  SI = [zeros(size(S(1:NE,1:NE))),sei*S(1:NE,NE+1:NN);zeros(size(S(NE+1:NN,1:NE))),sii*S(NE+1:NN,NE+1:NN)]; 
  
  RefTime = zeros(size(V_0s));
  SpikeRcd = [];
  oSp = sparse(size(V_0s,1),size(V_0s,2));
  
  for tInd = 1:steps
      %% nan represnets ref
      CurrentV = Vs(:,tInd);
      CurrentI_E = I_Es(:,tInd);
      CurrentI_I = I_Is(:,tInd);
      RefTime(isnan(CurrentV)) = RefTime(isnan(CurrentV)) + dt;
      CurrentV(RefTime>=t_ref) = 0;
      RefTime(RefTime>=t_ref) = 0;
%    dV = (dt*(-Vs(:,ii)/g_L) +(fs.*R(:,ii)+SE*spike(:,ii)-SI*spike(:,ii))/g_L ).*Ref(:,ii);
      dV = dt*(-CurrentV/g_L + CurrentI_E - CurrentI_I) ;
      oI_E = (CurrentI_E + (fs.*R(:,tInd) + SE*oSp)/tau_E) *exp(-dt/tau_E);
      oI_I = (CurrentI_I +                  SI*oSp /tau_I) *exp(-dt/tau_I);

      oV = CurrentV + dV; % the rough output
      oSp = sparse(double(oV>=1)); 
      
      oV(oV>=1) = nan;
      Vs(:,tInd+1) = oV;
      I_Es(:,tInd+1) = oI_E;
      I_Is(:,tInd+1) = oI_I;
      SpikeRcd = [SpikeRcd;[find(oSp),ones(size(find(oSp)))*tInd]] ; 
  end
  %% statistics
  spike = zeros(size(Vs));
  for spInd = 1:size(SpikeRcd,1)
      spike(SpikeRcd(spInd,1),SpikeRcd(spInd,2)) = 1;
  end
  % use norm_rate function to calculate firing rate m(t) for the whole network 
  pieces = floor(T/500); %500ms is the time length for each piece.
  %As in norm_rate function each entry in "rates" cover consecutive 2 pieces, actually the unit
  %for rates is # of spikes for 600 neurons per second.
  %rates = norm_rate(t,spike,pieces);
  SpikeCount = SpikeCount_slide(t,spike,5,1);
  
  V_Nm = norm_voltage(t,Vs,V_0s,pieces);  
  SynIndex = norm_syn(t,spike,N);
end