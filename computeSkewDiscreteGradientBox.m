function C = computeSkewDiscreteGradientBox(dim, p)
if isscalar(p)
  p = p * ones(dim,1);
end
n = p+1;
N = prod(n);
C = cell(1,dim);
for k = 1:dim
  c = computeSkewDiscreteGradient1D(p(k));
  C{k} = kron(kron(speye(prod(n(dim:-1:k+1))), c), speye(prod(n(1:k-1)))) * n(k) / prod(n);
end
end
