function [Vs,GEs,GIs,spike,t,SpikeCount] = MC_Vseries_for_network_ninput_G_based_ff(dt,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,I_gate,S,se,si)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  tau_L = 20;
  tau_f = 1;
  tau_E = 1;
  tau_I = 4;
%   sample_proportion = 1/5;
 
  N = length(V_0s); % number of neurons
    [A,~,Times] = find(eventsN);
    R = sparse(A,floor(Times/dt)+1,ones(size(A)),N,length(t));
    
  t_ref = 2;
 
  Vs = zeros(N,length(t));
  Vs(:,1) = V_0s;

  GEs = zeros(size(Vs)); GIs = zeros(size(Vs)); Gfs = zeros(size(Vs));
  GEs(:,1) = GE_0s;         GIs(:,1) = GI_0s;         Gfs(:,1) = Gf_0s;
  
  [row1, column1, weight1] = find(S>0);
  [row2, column2, weight2] = find(S<0);
  NN = length(V_0s);
  SE = se*sparse(row1, column1,weight1,NN,NN);
  SI = si*sparse(row2, column2,weight2,NN,NN); 
  
  VE = 14/3;
  VI = -2/3;
  
  RefTime = zeros(size(V_0s));
  SpikeRcd = [];
  oSp = sparse(size(V_0s,1),size(V_0s,2));
  for tInd = 1:steps
      %% nan represnets ref
      CurrentV = Vs(:,tInd);
      RefTime(isnan(CurrentV)) = RefTime(isnan(CurrentV)) + dt;
      CurrentV(RefTime>=t_ref) = 0;
      RefTime(RefTime>=t_ref) = 0;
      
   dV  = dt*(-CurrentV/tau_L - (Gfs(:,tInd)+GEs(:,tInd)).*(CurrentV-VE) - GIs(:,tInd).*(CurrentV-VI) - fs*I_gate);
   dGf = dt*1/tau_f*(-Gfs(:,tInd)) + fs.*R(:,tInd)/tau_f; % wrong for f/tau
   dGE = dt*1/tau_E*(-GEs(:,tInd)) + SE * oSp/tau_E;
   dGI = dt*1/tau_I*(-GEs(:,tInd)) + SI * oSp/tau_I;
   oV = CurrentV + dV; % the rough output
   oSp = sparse(double(oV>=1));
   
   oV(oV>=1) = nan;
   Vs(:,tInd+1) = oV;
   Gfs(:,tInd+1) = Gfs(:,tInd) + dGf;
   GEs(:,tInd+1) = GEs(:,tInd) + dGE;
   GIs(:,tInd+1) = GIs(:,tInd) + dGI;
   SpikeRcd = [SpikeRcd;[find(oSp),ones(size(find(oSp)))*tInd]] ; 
  end
%% statistics
  spike = zeros(size(Vs));
  for spInd = 1:size(SpikeRcd,1)
      spike(SpikeRcd(spInd,1),SpikeRcd(spInd,2)) = 1;
  end
  % use norm_rate function to calculate firing rate m(t) for the whole network 
  %pieces = floor(T/20); %500ms is the time length for each piece.
  %As in norm_rate function each entry in "rates" cover consecutive 2 pieces, actually the unit
  %for rates is # of spikes for 600 neurons per second.
  SpikeCount = SpikeCount_slide(t,spike,5,1); %norm_rate(t,spike,pieces);
end