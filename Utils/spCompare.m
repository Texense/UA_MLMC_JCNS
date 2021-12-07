%% spike comparison of two 1D spike-time series
% Assume spike2 is the more accurate one, treat it as template.
% spike1,2 should be both row vectors
% ZCX 03/11/21
function SpDiff = spCompare(spike12, spike22)
l1 = length(spike12);
l2 = length(spike22);
if l1 == l2 
    SpDiff = spike12 - spike22;
elseif max(l1,l2)>10000 % too large, just use simple solutions
    Sp1Short = l1(1:min(l1,l2));
    Sp2Short = l2(1:min(l1,l2));
    SpDiff = Sp1Short - Sp2Short;
else
    SpDiffMat = repmat(spike12',1,l2) - repmat(spike22,l1,1);
    [~,CompInd] = min(abs(SpDiffMat));
    SpDiff = zeros(size(spike22));
    for SpInd = 1:l2
        SpDiff(SpInd) = SpDiffMat(CompInd(SpInd),SpInd);
    end
end

end
