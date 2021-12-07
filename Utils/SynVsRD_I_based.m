function [V_diff,RelaDiff,Syn] = SynVsRD_I_based(T,ss1,ss2)
%T = 5000*10;

% N = 256;NE = 128;NI= 128;
N = 600;NE = 300;NI= 300;
V_0s = rand(N,1);
lambda = 0.55;
tau = 5; 

% I_0 = lambda*f*tau;
fs = [0.07*ones(NE,1);0.07*ones(NI,1)];
% see = 0.007;sei = 0.007;
% sie=0.007;   sii=0.007;

% pp = 1;
% pE = pp; pI = pp;
% SE = 1*((rand(N,NE))<pE);
% SI = 1*((rand(N,NI))<pI);
% S = [SE,SI];
S = ones(N);

for jj = 1:N
    S(jj,jj)=0;
end

Nevent = floor(T*lambda*2);
eventsN = zeros(N,Nevent);
for ii = 1:N
event = Poisson_Process_2(lambda,T);
eventsN(ii,1:length(event))=event;
end
clear ii
dt = 1*ones(1,7);
for ii = 1:length(dt)
   dt(ii) = dt(1)*2^(-ii+1); 
end

rate_ns = [];
V_ns = [];


tic
[Vs6,spike6,t6,rate_ns(:,7),V_ns(:,:,7),Syn6] = ...
    MC_Vseries_for_network_ninput_I_based_norms_Syn(dt(5),T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI);
toc
tic
[Vs7,spike7,t7,rate_ns(:,8),V_ns(:,:,8),Syn7] = ...
    MC_Vseries_for_network_ninput_I_based_norms_Syn(dt(6),T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI);
toc

V_diff = abs(Vs6(:,length(t6)-100) - Vs7(:,length(t6)-100));

T_sam = 500;
rates1 = norm_rate(t6,spike6,floor(T/T_sam));
rates_drop1 = rates1(1:end-1)*(1000/T_sam)/N;
rates2 = norm_rate(t7,spike7,floor(T/T_sam));
rates_drop2 = rates2(1:end-1)*(1000/T_sam)/N;
RelaDiff = abs(mean(rates_drop1)-mean(rates_drop2))/abs(mean(rates_drop1)+mean(rates_drop2));


Syn = (Syn6+Syn7)/2;

end




