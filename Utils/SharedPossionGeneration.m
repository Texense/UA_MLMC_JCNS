%% Generate N Possion Spike trains with shared input

function eventsN = SharedPossionGeneration(lambda1,ShareFactor,T,N)

%lambda1 = 0.55;
%ShareFactor = 0.0;
Nevent    = floor(T*lambda1*2);
NSHAEvent = floor(ShareFactor*Nevent); % number of shared input 
%NIIDEvent = Nevent - NSHAEvent;        % number of iid events

eventsN = zeros(N,Nevent);
SharedEvent = Poisson_Process_2(lambda1*ShareFactor,T);
for ii = 1:N
IIDevent = Poisson_Process_2(lambda1*(1-ShareFactor),T);
FinalEvent = sort([SharedEvent;IIDevent]);
eventsN(ii,1:length(FinalEvent)) = FinalEvent';
end
eventsN(:,1) = []; % get rid of the first column, all zero.
end