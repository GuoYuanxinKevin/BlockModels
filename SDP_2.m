function [Label, Agr] = SDP_2(A,n)
  %Label is the recovered community labels
  %Agr is the agreement between the recovered labels
  %and the true community labels.
  %0 <= Agr <= 1, where Agr == 1 implies exact recovery.
  bt = -[ones(1,n) 0]';
  B = A .* (A * A);
  ct = vec(-B);
  %B is the triangle matrix, i.e., 
  %B_ij = the number of triangles containing edge (i,j)
  At = [zeros(n^2,n) ones(n^2,1)];
  At(1:n^2+n+1:n^3) = 1;
  At = -At;
  K.s = n;
  [x,~,~] = sedumi(At,bt,ct,K);
  %Use SeDuMi to solve the (dual problem of the) SDP
  %Recall the SDP aims to recover the outer product of 
  %the true community label with itself.
  x = reshape(x,n,n);
  [Label,~] = eigs(x,1);
  Label = normr(Label);
  %Assign communities to the nodes by the sign of
  %the corresponding entry in the eigenvector 
  %associated to the largest eigenvalue.
  Agr = abs([ones(1,n/2) -1*ones(1,n/2)] * Label / n);
end