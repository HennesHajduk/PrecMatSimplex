function C = computeSkewDiscreteGradientSimplex(dim, p)
N = nchoosek(p+dim,p);
C = cell(1,dim);
alpha = getMultiindices(dim,p);
s = 0.5;
for k = 1:dim-1
  s = s / (p+k);
end
t = s / dim;
for i = 1:N
  ai = alpha(i,:);
  if ai(1) > 0
    aj = ones(dim,1) * ai + [-ones(dim,1), eye(dim)];
    for k = 1:dim
      j = getNodeFromMultiindex(dim, p, aj(k,:));
      for l = 1:dim
        C{l}(i,j) = t;
      end
    end
  else
    for l = 1:dim
      C{l}(i,i) = s;
    end
  end
  for k = 1:dim
    if ai(k+1) > 0
      j = getNodeFromMultiindex(dim, p, ai - [-1, k == (1:dim)]);
      for l = 1:dim
        C{l}(i,j) = -t;
      end
      for l = 1:dim
        if l == k, continue; end
        j = getNodeFromMultiindex(dim, p, ai + [0, l == (1:dim)] - [0, k == (1:dim)]);
        C{k}(i,j) = -t;
        C{l}(i,j) = t;
      end
    else
      C{k}(i,i) = C{k}(i,i) - s;
    end
  end
end
c = zeros(N,dim);
for k = 1:dim
  for j = 1:N
    if alpha(j,1) == 0
      c(j,k) = 2*s;
    end
    if alpha(j,k+1) == 0
      c(j,k) = c(j,k) - 2*s;
    end
  end
end
for k = 1:dim
  assert( max(abs(sum(C{k},1) - c(:,k)')) < 1e-15, [ 'Column sum wrong' ...
          ' on simplex for p = ' num2str(p) ', dim = ' num2str(dim) ] );
end
end
