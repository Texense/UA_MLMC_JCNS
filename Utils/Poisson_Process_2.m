function arrt = Poisson_Process_2(lambda,T)
npoints = poissrnd(lambda*T);
% When event number N is fixed, the event time are distributed averagely.
if (npoints>0)
  arrt = [0; sort(rand(npoints,1)*T)];
else
  arrt = [];
end
end