% inh_MC_Vend_for_single_neuron.m

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

function V_end = inh_MC_Vend_for_single_neuron(dt,T,V_0,I_0,event,f,I_inh)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms, we evolve 1s.
  g_L = 20;
  
  tau = 5;
  
  
%   sample_proportion = 1/5;
  R = zeros(size(t));
  V = zeros(size(t));
  V(1) = V_0;
  I = zeros(size(t));
  I(1) = I_0;
 

  for jj=1:length(event)
     R(floor(event(jj)/dt)+1) = R(floor(event(jj)/dt)+1)+1;   
  end

  spike = [];
 
  for ii = 1:steps
   
   dV = dt*(-V(ii)/g_L + I(ii) - I_inh);
   dI = dt*1/tau*(-I(ii)) + f*R(ii); % wrong for f/tau
   VV = V(ii) + dV;
     if VV > 1
       spike = [spike, ii*dt];
       V(ii+1) = VV - floor(VV);
     else V(ii+1) = VV;
     end
   
   I(ii+1) = I(ii) + dI;
  end
  V_end = V(steps+1);
  
end


% figure(1)
% plot(t,V)
% % figure(2)
% % plot(t,I)
