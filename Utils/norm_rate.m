function rates = norm_rate(t,spike,pieces)
  N = floor(length(t)/pieces);
  rates = zeros(1,pieces);
  for count = 1:pieces-1
      rates(count) = sum(sum(spike(:, 1+(count-1)*N:1+(count+1)*N)));
  end
end