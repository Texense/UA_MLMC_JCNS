%% JCNS Fig2
%% Recurrent Network; Ring network; Spike added to V
% You want to check if Same or different input are used for the samples
%% Setups
CurrentFolder = pwd;
addpath([CurrentFolder '/Utils'])

T = 1000;
%dt1 = 0.05; dt2 = dt1/2;

N = 600;NE = 300;NI= 300;
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
NeighborScale = 10; % 60 30 10 5
%Neighbor 5  10 30 60
N_conn = floor(NE/NeighborScale);
S_SWconn = zeros(NE);
for postInd = 1:NE
    for preInd = 1:NE
        if mod(postInd-preInd,NE)<=N_conn/2 | mod(postInd-preInd,NE)>= NE-N_conn/2
            S_SWconn(postInd,preInd) = 1;
        end
    end
end
S = repmat(S_SWconn,2,2);

for jj = 1:N
    S(jj,jj)=0;
end

ss1 = 0.003;%*NE/N_conn;
ss2 = 0.003;%*NE/(N_conn);


N_Sample = 32;
dtAll = [0.1,0.05]; Ndt = 2;
SharedFacNum = 11;
SharedProp = linspace(0,1,SharedFacNum);
%% Cluster Def
% creat a 10-hr parallel
cluster = gcp('nocreate');
if isempty(cluster)
    cluster = parpool([4 64]);
    cluster.IdleTimeout = 1200;
end

%% 10 SharedProp; Each include 50 pairs of coupled simulations 

% Vs    = cell(N_Sample,SharedFacNum,2);
% spike = cell(N_Sample,SharedFacNum,2);
% SpikeCount = cell(N_Sample,SharedFacNum,2);
% t     = cell(N_Sample,SharedFacNum,2);
% spiketime  = cell(N_Sample,SharedFacNum,2);

for ShareInd = 1:SharedFacNum
    ShareFactor = SharedProp(ShareInd);
    Nevent    = floor(T*lambda*2);
    NSHAEvent = floor(ShareFactor*Nevent); % number of shared input
    NIIDEvent = Nevent - NSHAEvent;        % number of iid events
    
    tic
    VsAll    = cell(N_Sample,2);
    spikeAll = cell(N_Sample,2);
    SpikeCountAll = cell(N_Sample,2);
    tAll     = cell(N_Sample,2);
    spiketimeAll  = cell(N_Sample,2);
    
    % Make up input event by mixing shared input and iid input
    eventsN = zeros(N,Nevent);
    SharedEvent = Poisson_Process_2(lambda*ShareFactor,T);
    for ii = 1:N
        IIDevent = Poisson_Process_2(lambda*(1-ShareFactor),T);
        FinalEvent = sort([SharedEvent;IIDevent]);
        eventsN(ii,1:length(FinalEvent)) = FinalEvent';
    end
    
    parfor SamInd = 1:N_Sample
        V_0s = rand(N,1); % IC of each sample is different
        

        
        for dtInd = 1:Ndt
            dt = dtAll(dtInd);
            
            [VsAll{SamInd, dtInd},spikeAll{SamInd, dtInd},...
             tAll{SamInd, dtInd}, SpikeCountAll{SamInd, dtInd},~,~] = ...
                MC_Vseries_for_network_ninput_I_based_norms_Syn(dt,T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI);
            
            [Sp_n,Sp_t] = find(spikeAll{SamInd, dtInd} == 1);
            spiketimeAll{SamInd, dtInd} = (Sp_t-1)*dt;
        end
        
    end
    
    for SamInd = 1:N_Sample
        Vs    = VsAll(SamInd,:);
        spike = spikeAll(SamInd,:);
        SpikeCount = SpikeCountAll(SamInd,:);
        t     = tAll(SamInd,:);
        spiketime = spiketimeAll(SamInd,:);
        save([CurrentFolder sprintf('/HPCData/fig2JCNS_ShF%.1f_Sam%d_SameInp.mat',ShareFactor,SamInd)],'Vs', 'spike','SpikeCount','t','spiketime',...
                                                                                               'ShareFactor','SharedProp','SamInd','eventsN')
    end
    toc
end
%% Trajs
% 
% fig2JCNS = ws2struct();
% save([CurrentFolder '/HPCData/fig2JCNS.mat'],'fig2JCNS','-v7.3')