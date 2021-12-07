%%this code evolve dynamics for a 600-neuron recuurent network. We first
%%randomly choose a initial condition of Vs for each neuron, then generate
%%uncorrelated Possion input spike trains to each neuron.We use the same initial
%%condition and the same input to evolve the network dynamics on different
%%time grids:dt=1ms to dt=2^-6ms.Output results include spike trains and
%%firing rates.
function [spiketime0,spiketime1,spiketime2,spiketime3,spiketime4] = Recurrent_multilevel_A(V_0s,N)
%Parameters
T = 5000;%simulation time

%N = 600;
NE = floor(N/2);NI= N-NE;
%V_0s = rand(N,1);%randomly chosen I.C.

lambda = 0.55;
tau = 5; 

fs = [0.07*ones(NE,1);0.07*ones(NI,1)];
see = 0.009;sei = 0.007;% connectivity
sie=0.007;   sii=0.009;

% exclude self recurrent
S = ones(N);
for jj = 1:N
    S(jj,jj)=0;
end

%generate uncorrelated Poission input to each neuron
Nevent = floor(T*lambda*2);
eventsN = zeros(N,Nevent); % eventsN records inputs to each neuron
for ii = 1:N
    event = Poisson_Process_2(lambda,T);
    eventsN(ii,1:length(event))=event;
end
clear ii

%time grids
% dt = 1*ones(1,7);
% for ii = 1:length(dt)
%     dt(ii) = dt(1)*2^(-ii+1); 
% end
dt = [0.4 0.2 0.1 0.05 0.025];

rate_ns = [];
%% evolve the network dynamics on different level. get spiketime matrix, and rate vector
tic
[Vs0,spike0,t0,rate_ns(:,1)] = ...
    MC_Vseries_for_network_ninput_I_based_norms(dt(1),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
toc

tic
[Vs1,spike1,t1,rate_ns(:,2)] = ...
    MC_Vseries_for_network_ninput_I_based_norms(dt(2),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
toc
tic
[Vs2,spike2,t2,rate_ns(:,3)] = ...
    MC_Vseries_for_network_ninput_I_based_norms(dt(3),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
toc

tic
[Vs3,spike3,t3,rate_ns(:,4)] = ...
    MC_Vseries_for_network_ninput_I_based_norms(dt(4),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
toc

tic
[Vs4,spike4,t4,rate_ns(:,5)] = ...
    MC_Vseries_for_network_ninput_I_based_norms(dt(5),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
toc

% tic
% [Vs5,spike5,t5,rate_ns(:,6)] = ...
%     MC_Vseries_for_network_ninput_I_based_norms(dt(6),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
% toc
% tic
% [Vs6,spike6,t6,rate_ns(:,7)] = ...
%     MC_Vseries_for_network_ninput_I_based_norms(dt(7),T,V_0s,eventsN,fs,S,see,sei,sie,sii,NE,NI);
% toc
%% get specific spike time points for each level
[Sp_n0,Sp_t0] = find(spike0 == 1);
spiketime0 = (Sp_t0-1)*dt(1);
[Sp_n1,Sp_t1] = find(spike1 == 1);
spiketime1 = (Sp_t1-1)*dt(2);
[Sp_n2,Sp_t2] = find(spike2 == 1);
spiketime2 = (Sp_t2-1)*dt(3);
[Sp_n3,Sp_t3] = find(spike3 == 1);
spiketime3 = (Sp_t3-1)*dt(4);
[Sp_n4,Sp_t4] = find(spike4 == 1);
spiketime4 = (Sp_t4-1)*dt(5);
% [Sp_n5,Sp_t5] = find(spike5 == 1);
% spiketime5 = (Sp_t5-1)*dt(6);
% [Sp_n6,Sp_t6] = find(spike6 == 1);
% spiketime6 = (Sp_t6-1)*dt(7);
end
