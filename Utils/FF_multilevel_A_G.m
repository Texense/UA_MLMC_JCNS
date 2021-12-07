% FF_multilevel_A_G.m

% This file is part of UA_MLMC_JCNS, a collection of
% numerical experiments on multilevel Monte Carlo for
% spiking neuron networks.   It accompanies the
% paper "Multilevel Monte Carlo for Cortical Circuit Models,"
% to appear in the Journal of Computational Neuroscience.
% 
% Copyright (C) 2021 by Zhuo-Cheng Xiao
% <zx555@nyu.edu>
% 
% This program is free software; you can redistribute
% it and/or modify it under the terms of the GNU
% General Public License as published by the Free
% Software Foundation; either version 2 of the
% License, or (at your option) any later version.
% 
% This program is distributed in the hope that it
% will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General
% Public License along with this program; if not, write
% to the Free Software Foundation, Inc., 51 Franklin
% Street, Fifth Floor, Boston, MA 02110-1301 USA.

function [spiketime0,spiketime1,spiketime2,spiketime3,spiketime4] = FF_multilevel_A_G(V_0s,N)
T = 200;

%N = 10;
%V_0s = zeros(N,1);
Gf_0s = zeros(N,1);
GE_0s = zeros(N,1);
GI_0s = zeros(N,1);


lambda = 1;
tau = 5;
I_gate = lambda*tau;
% I_0 = lambda*f*tau;
fs = 0.05*ones(N,1);
se = 0.02;si=0.1;
S11 = [0,0,0,0,0;
         1,0,0,0,0;
         -1,1,0,0,0;
         -1,1,1,0,0;
         1,-1,1,-1,0];
S22 = [0,0,0,0,0;
         -1,0,0,0,0;
         1,-1,0,0,0;
         0,-1,1,0,0;
         1,0,1,-1,0];
S12 = [0,-1,1,0,1;
          0,1,1,0,-1;
         0,1,0,0,1;
         0,-1,1,0,1;
         0,1,0,1,0];
S = [S11,zeros(5);S12, S22];

Nevent = floor(T*lambda*2);
eventsN = zeros(N,Nevent);
for ii = 1:N
event = Poisson_Process_2(lambda,T);
eventsN(ii,1:length(event))=event;
end
clear ii
dt0= 0.2/2;
dt1 = 0.1/2;
dt2 = 0.05/2;
dt3 = 0.025/2;
dt4 = 0.0125/2; 



[Vs0,GE0s,GI0s,spike0,t0,rates0] = MC_Vseries_for_network_ninput_G_based_ff(dt0,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,0,S,se,si);

[Vs1,GE1s,GE1s,spike1,t1,rates1] = MC_Vseries_for_network_ninput_G_based_ff(dt1,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,0,S,se,si);

[Vs2,GE2s,GI2s,spike2,t2,rates2] = MC_Vseries_for_network_ninput_G_based_ff(dt2,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,0,S,se,si);

[Vs3,GE3s,GI3s,spike3,t3,rates3] = MC_Vseries_for_network_ninput_G_based_ff(dt3,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,0,S,se,si);

tic
[Vs4,GE4s,GI4s,spike4,t4,rates4] = MC_Vseries_for_network_ninput_G_based_ff(dt4,T,V_0s,Gf_0s,GE_0s,GI_0s,eventsN,fs,0,S,se,si);
toc


[Sp_n0,Sp_t0] = find(spike0 == 1);
spiketime0 = (Sp_t0-1)*dt0;
[Sp_n1,Sp_t1] = find(spike1 == 1);
spiketime1 = (Sp_t1-1)*dt1;
[Sp_n2,Sp_t2] = find(spike2 == 1);
spiketime2 = (Sp_t2-1)*dt2;
[Sp_n3,Sp_t3] = find(spike3 == 1);
spiketime3 = (Sp_t3-1)*dt3;
[Sp_n4,Sp_t4] = find(spike4 == 1);
spiketime4 = (Sp_t4-1)*dt4;
end

% % make graph
% [row, column,weight] = find(S');
% [row1, colunm1, weight1] = find(S'>0);
% [row2, colunm2, weight2] = find(S'<0);
% Graph = digraph(row, column, weight);
% figure(1);
% subplot(3,2,[1,3]);
% gg = plot(Graph);
% highlight(gg,row2, colunm2,'EdgeColor','b','LineWidth',1.5)
% highlight(gg,row1, colunm1,'EdgeColor','r','LineWidth',1.5)
% gg.XData = [0 1 2 3 4 4 5 5 6 7];
% gg.YData = [5 0 6 12 -1 10 0 8 2 9];
% axis([0 7 -2 13])
% text(7,12,'A');
% axis off
% subplot(3,2,5);
% 
% scatter(Sp_t3*dt3,Sp_n3,1,'r');
% axis([0 T 0 N])
% edges = -5:0.25:10; 
% subplot(3,2,2);histogram(spiketime1(1:length(spiketime0))-spiketime0,dt0*edges); 
% %axis([-2*dt0 6*dt0 0 60]);text(5.5*dt0,50,'B');
% subplot(3,2,4);histogram(spiketime2(1:length(spiketime1))-spiketime1,dt1*edges);
% %axis([-2*dt1 6*dt1 0 60]);text(5.5*dt1,50,'C');
% subplot(3,2,6);histogram(spiketime3(1:length(spiketime2))-spiketime2,dt2*edges);
% %axis([-2*dt2 6*dt2 0 60]);text(5.5*dt2,50,'D');
% xlabel('error of spike time(ms)')
% 
% %print(1,'-r600','-dpng','G_based_feedforward.png');
% %print(1,'-r600','-dpdf','G_based_feedforward.pdf');