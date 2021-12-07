function [Vs,Gfs,GEs,GIs,spike,t,rates] = MC_Vseries_for_network_ninput_G_based_norms(...
                                              dt1,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,I_gate,S,see,sei,sie,sii,NE,NI)
  steps = floor(T/dt1);
  t = 0:dt1:T; % time line. unit is ms.
  tau_L = 20;
  tau_f = 5;
  tau_E = 5;
  tau_I = 10;
%   tau_L = 20;
%   tau_f = 5;
%   tau_E = 50;
%   tau_I = 20;
  N = length(V_0s);
%   sample_proportion = 1/5;
  [A,B,Times] = find(eventsN);
  R = sparse(A,floor(Times/dt1)+1,ones(size(A)),N,length(t));

   % number of neurons
%   t_ref = 5;
  t_ref = 5;
  N_ref = floor(t_ref/dt1);
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
   
   dV = (dt1*(-Vs(:,ii)/tau_L - (Gfs(:,ii)+GEs(:,ii)).*(Vs(:,ii)-VE) - GIs(:,ii).*(Vs(:,ii)-VI) - fs*I_gate)).*Ref(:,ii);
   dGf = dt1*1/tau_f*(-Gfs(:,ii)) + fs.*R(:,ii)/tau_f; % wrong for f/tau
   dGE = dt1*1/tau_E*(-GEs(:,ii)) + SE*spike(:,ii)/tau_E;
   dGI = dt1*1/tau_I*(-GEs(:,ii)) + SI*spike(:,ii)/tau_I;
   VV = Vs(:,ii) + dV; % the rough output
   V_pls = mod(subplus(VV),1);
   spike(:,ii+1) = subplus(VV)-V_pls;
   Ref(:,ii+1:ii+N_ref) = Ref(:,ii+1:ii+N_ref) - repmat(spike(:,ii+1),[1,N_ref]);
   Vs(:,ii+1) = min(V_pls,VV);
   Gfs(:,ii+1) = Gfs(:,ii) + dGf;
   GEs(:,ii+1) = GEs(:,ii) + dGE;
   GIs(:,ii+1) = GIs(:,ii) + dGI;
  end
  pieces = 500;
  rates = norm_rate(t,spike,pieces);
%   V_Nm = norm_voltage(t,Vs,V_0s,pieces);
%   G_Nm = norm_conductance(t,Gfs,GEs,GIs,V_0s,pieces);
end