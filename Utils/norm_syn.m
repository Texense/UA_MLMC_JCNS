function SynIndex = norm_syn(t,spike,N)

SpCountTime = sparse(sum(spike)');
Win = 2.5;
dt0 = t(2)-t(1);
WinGrid = floor(Win/dt0);
SynKernel = spdiags(ones(length(t),2*WinGrid+1),-WinGrid:WinGrid,length(t),length(t));

SynIndex = SpCountTime'*SynKernel*SpCountTime / N / sum(SpCountTime);



end
