function C = computeSkewDiscreteGradientPrism(pTri, p1D)
n1D = p1D + 1;
nTri = nchoosek(pTri+2,pTri);
C = cell(1,3);
cTri = computeSkewDiscreteGradientSimplex(2, pTri);
C{1} = kron(speye(n1D), cTri{1}) / n1D;
C{2} = kron(speye(n1D), cTri{2}) / n1D;
C{3} = kron(computeSkewDiscreteGradient1D(p1D), speye(nTri)) / ((pTri+1)*(pTri+2));
s = 1/(pTri+1);
cTri = zeros(nTri,2);
a = getMultiindices(2, pTri);
for k = 1:2
  for j = 1:nTri
    if a(j,1) == 0
      cTri(j,k) = s;
    end
    if a(j,k+1) == 0
      cTri(j,k) = cTri(j,k) - s;
    end
  end
end
tol = 1e-15;
c = [ [kron(ones(n1D,1), cTri(:,1)), kron(ones(n1D,1), cTri(:,2))] / n1D, ...
      kron([-1 zeros(1, n1D-2) 1]', ones(nTri,1)) / ((pTri+1)*(pTri+2)) ];
for k = 1:3
  tmp = C{k};
  str = [' on prism for p = ' num2str([pTri p1D])];
  assert(max(abs(sum(tmp,2))) < tol, ['Rowsum wrong' str]);
  assert(max(abs(sum(tmp,1) - c(:,k)')) < tol, ['Column sum wrong' str]);
  tmp = tmp - diag(diag(tmp));
  assert(max(max(abs(tmp + tmp'))) < tol, ['No skew symmetry' str]);
end
end
