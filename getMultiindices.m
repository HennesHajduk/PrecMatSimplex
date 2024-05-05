function ret = getMultiindices(dim, p)
N = nchoosek(p+dim,p);
ret = zeros(N,dim+1);
a = [p zeros(1,dim)];
ctr = 1;
if dim == 1
  ret = [p:-1:0; 0:p]';
elseif dim == 2
  for j = 0:p
    for i = 0:p-j
      ret(ctr,:) = a;
      ctr = ctr+1;
      a(1:2) = a(1:2) + [-1 1];
    end
    a(2) = 0;
    a(3) = a(3)+1;
    a(1) = p - a(3);
  end
elseif dim == 3
  for k = 0:p
    for j = 0:p-k
      for i = 0:p-k-j
        ret(ctr,:) = a;
        ctr = ctr+1;
        a(1:2) = a(1:2) + [-1 1];
      end
      a(2) = 0;
      a(3) = a(3)+1;
      a(1) = p - sum(a(3:4));
    end
    a(2:3) = 0;
    a(4) = a(4)+1;
    a(1) = p - a(4);
  end
else
  error('Spatial dimension not supported.');
end
end
