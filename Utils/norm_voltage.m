function V_Nm = norm_voltage(t,Vs,V_0s,pieces)
  N = floor(length(t)/pieces);
  V_Nm = zeros(length(V_0s),pieces);
  for count = 1:pieces-1
      V_Nm(:,count) = Vs(:,N*count);
  end
end