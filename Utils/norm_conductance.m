function G_Nm = norm_conductance(t,Gfs,GEs,GIs,V_0s,pieces)
  N = floor(length(t)/pieces);
  G1 = zeros(length(V_0s),pieces);
  G2 = zeros(length(V_0s),pieces);
  G3 = zeros(length(V_0s),pieces);
  for count = 1:pieces-1
      G1(:,count) = Gfs(:,N*count);
      G2(:,count) = GEs(:,N*count);
      G3(:,count) = GIs(:,N*count);
  end
  G_Nm = [G1;G2;G3];
end