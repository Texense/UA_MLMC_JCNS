%% JCNS Fig2
%% Recurrent Network; Ring network; Spike added to V
% You want to check if Same or different input are used for the samples
%% Setups
CurrentFolder = pwd;
addpath([CurrentFolder '/Utils'])
HPCPath = [CurrentFolder '/HPCData'];

T = 1000;

% N = 256;NE = 128;NI= 128;
N = 600;NE = 300;NI= 300;


lambda = 0.55;

% I_0 = lambda*f*tau;
fs = [0.07*ones(NE,1);0.07*ones(NI,1)];
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

N_Sample = 32;
dtAll = [0.1, 0.05];

%% Cluster Def
% creat a 10-hr parallel
cluster = gcp('nocreate');
if isempty(cluster)
    cluster = parpool([4 64]);
    cluster.IdleTimeout = 1200;
end

%% Each (ss1,ss2) include 32 pairs of coupled sample (all)

for ss = 0.005:0.005:0.05
    see = ss;
    sei = ss*0.9;
    sie = ss*0.9;
    sii = ss;
    
    
    tic
    VsAll    = cell(N_Sample,2);
    spikeAll = cell(N_Sample,2);
    SpikeCountAll = cell(N_Sample,2);
    tAll     = cell(N_Sample,2);
    spiketimeAll  = cell(N_Sample,2);
    
    % Make up input event by mixing shared input and iid input
    Nevent = floor(T*lambda*2);
    eventsN = zeros(N,Nevent);
    for ii = 1:N
        event = Poisson_Process_2(lambda,T);
        eventsN(ii,1:length(event))=event;
    end
    
    WinSize  = 40;
    WinSlide = 20;
    parfor SamInd = 1:N_Sample
        V_0s = rand(N,1); % IC of each sample is different
        

        
        for dtInd = 1:Ndt
            dt = dtAll(dtInd);
            
            [VsAll{SamInd, dtInd},spikeAll{SamInd, dtInd},...
             tAll{SamInd, dtInd}, SpikeCountAll{SamInd, dtInd},~,~] = ...
                MC_Vseries_for_network_ninput_I_based_norms_Syn(dt,T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
            
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
        save([CurrentFolder sprintf('/HPCData/fig2JCNS_ss%.3f_Sam%d_SameInp.mat',ss,SamInd)],'Vs', 'spike','SpikeCount','t','spiketime',...
                                                                                               'ShareFactor','SharedProp','SamInd','eventsN')
    end
    toc
end

%% Trajs
% 
% fig2JCNS = ws2struct();
% save([CurrentFolder '/HPCData/fig2JCNS.mat'],'fig2JCNS','-v7.3')