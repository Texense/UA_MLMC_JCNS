function [rates_drop,RelaDiff,Syn] = RecurrentNetwork_couple(T,ss1,ss2,dt1,dt2,T_sam)
%T = 5000*10;

% N = 256;NE = 128;NI= 128;
N = 600;NE = 300;NI= 300;
V_0s = rand(N,1);
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
% event = Poisson_Process_2(lambda,T);
% for ii = 1:N
% eventsN(ii,1:length(event))=event;
% end


[~,spike1,t1,~,~,Syn1] = ...
    MC_Vseries_for_network_ninput_I_based_norms_Syn(dt1,T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI);

[~,spike2,t2,~,~,Syn2] = ...
    MC_Vseries_for_network_ninput_I_based_norms_Syn(dt2,T,V_0s,eventsN,fs,S,ss1,ss2,ss2,ss1,NE,NI);



%T_sam = 200;
rates1E = norm_rate(t1,spike1(1:NE,:),floor(T/T_sam));   rates_drop1E = rates1E(1:end-1)*(1000/T_sam)/NE;
rates1I = norm_rate(t1,spike1(NE+1:N,:),floor(T/T_sam)); rates_drop1I = rates1I(1:end-1)*(1000/T_sam)/NI;
rates_drop1 = [rates_drop1E;rates_drop1I];

rates2E = norm_rate(t2,spike2(1:NE,:),floor(T/T_sam));   rates_drop2E = rates2E(1:end-1)*(1000/T_sam)/NE;
rates2I = norm_rate(t2,spike2(NE+1:N,:),floor(T/T_sam)); rates_drop2I = rates2I(1:end-1)*(1000/T_sam)/NI;
rates_drop2 = [rates_drop2E;rates_drop2I];


rates_drop(:,:,1) = rates_drop1;
rates_drop(:,:,2) = rates_drop2;

RelaDiff = mean((rates_drop1-rates_drop2)./(rates_drop1+rates_drop2),'all');


Syn = (Syn1+Syn2)/2;

end