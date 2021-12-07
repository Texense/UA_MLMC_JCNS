function [Vs,Is,spike,t] = MC_Vseries_for_network_oneinput(dt,T,V_0s,I_0s,event,fs,I_gate,S)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  g_L = 20;
  tau = 5;
  
%   sample_proportion = 1/5;
  R = zeros(size(t));
  for jj=1:length(event)
     R(floor(event(jj)/dt)+1) = R(floor(event(jj)/dt)+1)+1;   
  end
  
  N = length(V_0s); % number of neurons
  spike = zeros(N,length(t));
 
  Vs = zeros(N,length(t));
  Vs(:,1) = V_0s;
  Is = zeros(N,length(t));
  Is(:,1) = I_0s;
  
  for ii = 1:steps
   
   dV = dt*(-Vs(:,ii)/g_L + Is(:,ii) - fs*I_gate);
   dI = dt*1/tau*(-Is(:,ii)) + (fs*R(ii) + S*spike(:,ii))/tau;
   VV = Vs(:,ii) + dV; % the rough output
   V_pls = mod(subplus(VV),1);
   spike(:,ii+1) = subplus(VV)-V_pls;
   Vs(:,ii+1) = min(V_pls,VV);
   Is(:,ii+1) = Is(:,ii) + dI;
  end
    
end