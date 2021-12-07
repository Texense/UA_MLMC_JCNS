%% main function to evolve recurrent network
function [Vs,spike,t,rates,V_Nm,SynIndex] = MC_Vseries_for_network_ninput_I_based_norms_Syn(dt,T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI,varargin)
  steps = floor(T/dt);
  t = 0:dt:T; % time line. unit is ms.
  g_L = 20;
  
  if nargin <=12
      WinSize  = 10;
      WinSlide = 5;
  else
      WinSize = varargin{1};
      WinSlide = varargin{2};
  end 
  

  N = length(V_0s);
% Read input spike time points in eventsN, and put the points in numerical
% iteration process
  [A,~,Times] = find(eventsN);
  R = sparse(A,floor(Times/dt)+1,ones(size(A)),N,length(t));
  %R = full(R);
  % add refractory period
  t_ref = 1;
 
  % matrix Vs records V for every neuron on every time step
  Vs = zeros(N,length(t));
  Vs(:,1) = V_0s;
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
      RefTime(isnan(CurrentV)) = RefTime(isnan(CurrentV)) + dt;
      CurrentV(RefTime>=t_ref) = 0;
      RefTime(RefTime>=t_ref) = 0;
%    dV = (dt*(-Vs(:,ii)/g_L) +(fs.*R(:,ii)+SE*spike(:,ii)-SI*spike(:,ii))/g_L ).*Ref(:,ii);
      dV = dt*(-CurrentV/g_L) +(fs.*R(:,tInd)+SE*oSp-SI*oSp);
      oV = CurrentV + dV; % the rough output
      oSp = sparse(double(oV>=1)); 
      
      oV(oV>=1) = nan;
      Vs(:,tInd+1) = oV;
      
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
  rates = SpikeCount_slide(t,spike,WinSize,WinSlide);
  
  V_Nm = norm_voltage(t,Vs,V_0s,pieces);  
  SynIndex = norm_syn(t,spike,N);
end