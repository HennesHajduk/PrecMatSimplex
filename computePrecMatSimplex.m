function C = computePrecMatSimplex(dim, p)
N = nchoosek(p+dim,p);
C = cell(1,dim);
multiindices = getMultiindices(dim,p);
volRefSimplex = 1 / factorial(dim);
M_L = volRefSimplex / N;
gradBary = [-ones(1,dim); eye(dim)];
for m = 1:dim
  C{m} = sparse(N,N);
  for j = 1:N
    alpha = multiindices(j,:);
    for k = 1:dim+1
      for l = 1:dim+1
        beta = alpha;
        beta(k) = beta(k) - 1;
        beta(l) = beta(l) + 1;
        if (beta(l) > p || beta(k) < 0)
          continue;
        end
        i = getNodeFromMultiindex(dim, p, beta);
        C{m}(i,j) = C{m}(i,j) + gradBary(k,m)  * (alpha(l) - (l==k) + 1);
      end
    end
  end
  C{m} = M_L * C{m};
end
end
