function testDiscreteGradients
%% Tensor-product elements with isotropic polynomial degrees
P = [56 31 15 7 3];
for dim = 1:5
  for p = 1:P(dim)
    C1 = computePrecMatBox(dim, p);
    C2 = computeSkewDiscreteGradientBox(dim, p);
    str = [' on tensor-product element for p = ' num2str(p) ', dim = ' num2str(dim)];
    runTests(C1, C2, str);
  end
end
fprintf('Isotropic tensor-product element tests ran without failures.\n');
%% Tensor-product elements with anisotropic polynomial degrees
P = [10, 5, 4, 3];
for dim = 2:5
  p = ones(1,dim);
  for i = 1:P(dim-1)^dim
    C1 = computePrecMatBox(dim, p);
    C2 = computeSkewDiscreteGradientBox(dim, p);
    str = [' on tensor-product element for p = ' num2str(p) ', dim = ' num2str(dim)];
    runTests(C1, C2, str);
    p(1) = p(1) + 1;
    for l = 1:dim
      if p(l) > P(dim-1)
        p(l) = 1;
        if l < dim
          p(l+1) = p(l+1) + 1;
        end
      end
    end
  end
end
fprintf('Anisotropic tensor-product element tests ran without failures.\n');
%% Isotropic simplices (anisotropy not allowed for Bernstein elements on simplices)
%% An additional test is performed within the called routine
for dim = 1:3 % if dim == 3, computation requires signifcant time for high polynomial degrees
  for p = 1:31
    C1 = computePrecMatSimplex(dim, p);
    C2 = computeSkewDiscreteGradientSimplex(dim, p);
    if dim == 1
      C3 = computeSkewDiscreteGradientBox(1, p);
      assert(isequal(C2{1}, C3{1}), 'Matrices on boxes and simplices must be the same in 1D.');
    end
    str = [' on simplex element for p = ' num2str(p) ', dim = ' num2str(dim)];
    runTests(C1, C2, str);
  end
end
fprintf('Simplex element tests ran without failures.\n');
%% Anisotropic prisms (triangles in (x,y)-plane times [0,1] in z-direction)
%% Tests are performed within the called routine
for p1D = 1:31
  for pTri = 1:31
    computeSkewDiscreteGradientPrism(pTri, p1D);
  end
end
fprintf('Anisotropic prism element tests ran without failures.\n');
end

function runTests(C1, C2, str)
tol = 1e-13;
dim = size(C1,3);
for k = 1:dim
  tmp = C2{k} + C2{k}';
  tmp = tmp - diag(diag(tmp));
  assert(max(max(abs(tmp))) < tol, ['No skew symmetry' str]);
  assert(max(abs(sum(C1{k}, 2))) < tol, ['Dense matrix rowsum wrong' str]);
  assert(max(abs(sum(C2{k}, 2))) < tol, ['Sparse matrix rowsum wrong' str]);
  assert(max(max(abs(sum(C1{k} - C2{k}, 1)))) < tol, ['Column sums wrong' str]);
end
end
