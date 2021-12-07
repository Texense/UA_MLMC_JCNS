function rates = SpikeCount_slide(t,spike,WinSize,slideSize)
  pieces = floor((t(end)-WinSize)/slideSize);
  NWinBin = floor(WinSize/(t(2)-t(1)));
  SlideBin = floor(slideSize/(t(2)-t(1)));
  rates = zeros(size(spike,1),pieces);
  for WinInd = 1:pieces
      rates(:,WinInd) = sum(spike(:, 1+(WinInd-1)*SlideBin : (WinInd-1)*SlideBin + NWinBin),2);
  end
end