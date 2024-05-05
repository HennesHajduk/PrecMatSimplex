function C = computePrecMatBox(dim, p)
if isscalar(p)
  p = p*ones(1,dim);
end
N = prod(p+1);
C = cell(1,dim);
for k = 1:dim
  c1D = computePrecMat1D(p(k));
  C{k} = kron( speye(prod(p(k+1:end)+1)), kron( c1D{1}, speye(prod(p(1:k-1)+1)) ) ) * (p(k)+1) / N;
end
end
