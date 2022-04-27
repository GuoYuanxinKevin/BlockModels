function [Loc, Adj] = Graph_Gen(n,p,q,d)
  %Assume by default the first n/2 nodes are in first community
  %and the last n/2 are in the second
  %p and q are the in-cluster and across-cluster probabilties
  %Assume p,q < 0.5
  mu = zeros(1,d+1);
  sig = eye(d+1);
  Loc = mvnrnd(mu,sig,n);
  Loc = normr(Loc);
  %Generate the locations of the n nodes
  %Nodes are uniformly dist'd on the d-dim'l sphere
  %Project multivariate Gaussian to the sphere
  Thres = [p*ones(n/2),q*ones(n/2); q*ones(n/2), p*ones(n/2)] - eye(n);
  CosTheta = max(Loc * Loc', 0);
  Adj = 0.5 * betainc(1-min(1,CosTheta.^2),d/2,0.5);
  Adj(Adj > Thres) = 0;
  Adj(Adj > 0) = 1;
  %The ratio of the surface area of a hyperspherical cap
  %to the entire surface area is 0.5 * betainc(sin(theta)^2,d/2,0.5)
  %where theta is the colatitude.
end
