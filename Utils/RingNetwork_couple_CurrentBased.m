function [rates_drop,RelaDiff,Syn] = RingNetwork_couple_CurrentBased(T,ss1,ss2,dt1,dt2,T_sam,ShareFactor,NeighborScale)
%T = 5000*10;

% N = 256;NE = 128;NI= 128;
N = 600;NE = 300;NI= 300;
V_0s = rand(N,1);
I_E0 = rand(N,1)*0.05;
I_I0 = rand(N,1)*0.05;
lambda = 0.55;
%tau = 5; 

% I_0 = lambda*f*tau;
fs = [0.07*ones(NE,1);0.07*ones(NI,1)];
% see = 0.007;sei = 0.007;
% sie=0.007;   sii=0.007;

% pp = 1;
% pE = pp; pI = pp;
% SE = 1*((rand(N,NE))<pE);
% SI = 1*((rand(N,NI))<pI);
% S = [SE,SI];
%NeighborScale = 30; % 50 for ring model
N_conn = floor(NE/NeighborScale);
S_SWconn = zeros(NE);
for postInd = 1:NE
    for preInd = 1:NE
        if mod(postInd-preInd,NE)<=N_conn/2 || mod(postInd-preInd,NE)>= NE-N_conn/2
            S_SWconn(postInd,preInd) = 1;
        end  
    end 
end
S = repmat(S_SWconn,2,2);

for jj = 1:N
    S(jj,jj)=0;
end

%ShareFactor = 0.0;
Nevent    = floor(T*lambda*2);
% NSHAEvent = floor(ShareFactor*Nevent); % number of shared input 
% NIIDEvent = Nevent - NSHAEvent;        % number of iid events

eventsN = zeros(N,Nevent);
SharedEvent = Poisson_Process_2(lambda*ShareFactor,T);
for ii = 1:N
IIDevent = Poisson_Process_2(lambda*(1-ShareFactor),T);
FinalEvent = sort([SharedEvent;IIDevent]);
eventsN(ii,1:length(FinalEvent)) = FinalEvent';
end


[~,~,~,spike1,t1,~,~,Syn1] = ...
    MLMC_LIF_Current_Based_Syn(dt1,T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI,...
                               I_E0,I_I0);


[~,~,~,spike2,t2,~,~,Syn2] = ...
    MLMC_LIF_Current_Based_Syn(dt2,T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI,...
                               I_E0,I_I0);



%T_sam = 200;
rates1E = norm_rate(t1,spike1(1:NE,:),floor(T/T_sam));   rates_drop1E = rates1E(1:end-1)*(1000/T_sam)/NE/2;
rates1I = norm_rate(t1,spike1(NE+1:N,:),floor(T/T_sam)); rates_drop1I = rates1I(1:end-1)*(1000/T_sam)/NI/2;
rates_drop1 = [rates_drop1E;rates_drop1I];

rates2E = norm_rate(t2,spike2(1:NE,:),floor(T/T_sam));   rates_drop2E = rates2E(1:end-1)*(1000/T_sam)/NE/2;
rates2I = norm_rate(t2,spike2(NE+1:N,:),floor(T/T_sam)); rates_drop2I = rates2I(1:end-1)*(1000/T_sam)/NI/2;
rates_drop2 = [rates_drop2E;rates_drop2I];


rates_drop(:,:,1) = rates_drop1;
rates_drop(:,:,2) = rates_drop2;

RelaDiff = mean((rates_drop1-rates_drop2)./(rates_drop1+rates_drop2),'all');


Syn = (Syn1+Syn2)/2;

end