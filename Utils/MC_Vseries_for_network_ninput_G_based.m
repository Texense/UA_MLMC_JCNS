function [Vs,Gfs,GEs,GIs,spike,t] = MC_Vseries_for_network_ninput_G_based(dt,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,I_gate,S,see,sei,sie,sii,NE,NI)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  tau_L = 20;
  tau_f = 1;
  tau_E = 3;
  tau_I = 2;
  N = length(V_0s);
%   sample_proportion = 1/5;
  [A,B,Times] = find(eventsN);
  R = sparse(A,floor(Times/dt)+1,ones(size(A)),N,length(t));

   % number of neurons
  t_ref = 5;
  N_ref = floor(t_ref/dt);
  spike = zeros(N,length(t));
    
  Vs = zeros(N,length(t)+N_ref);
  Vs(:,1) = V_0s;
  Ref = ones(size(Vs)); 
  GEs = zeros(size(Vs)); GIs = zeros(size(Vs)); Gfs = zeros(size(Vs));
  GEs(:,1) = GE_0s;         GIs(:,1) = GI_0s;         Gfs(:,1) = Gf_0s;
   % portion of Gfs in excitory
  NN = length(V_0s);
  SE = [see*S(1:NE,1:NE),zeros(size(S(1:NE,NE+1:NN)));sie*S(NE+1:NN,1:NE),zeros(size(S(NE+1:NN,NE+1:NN)))];
  SI = [zeros(size(S(1:NE,1:NE))),sei*S(1:NE,NE+1:NN);zeros(size(S(NE+1:NN,1:NE))),sii*S(NE+1:NN,NE+1:NN)]; 
  
  VE = 14/3;
  VI = -2/3;
  for ii = 1:steps
   
   dV = (dt*(-Vs(:,ii)/tau_L - (Gfs(:,ii)+GEs(:,ii)).*(Vs(:,ii)-VE) - GIs(:,ii).*(Vs(:,ii)-VI) - fs*I_gate)).*Ref(:,ii);
   dGf = dt*1/tau_f*(-Gfs(:,ii)) + fs.*R(:,ii); % wrong for f/tau
   dGE = dt*1/tau_E*(-GEs(:,ii)) + SE*spike(:,ii);
   dGI = dt*1/tau_I*(-GEs(:,ii)) + SI*spike(:,ii);
   VV = Vs(:,ii) + dV; % the rough output
   V_pls = mod(subplus(VV),1);
   spike(:,ii+1) = subplus(VV)-V_pls;
   Ref(:,ii+1:ii+N_ref) = Ref(:,ii+1:ii+N_ref) - repmat(spike(:,ii+1),[1,N_ref]);
   Vs(:,ii+1) = min(V_pls,VV);
   Gfs(:,ii+1) = Gfs(:,ii) + dGf;
   GEs(:,ii+1) = GEs(:,ii) + dGE;
   GIs(:,ii+1) = GIs(:,ii) + dGI;
  end
    
end