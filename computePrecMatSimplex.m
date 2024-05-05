% For the reference simplex $\hat K$ in dim dimensions this routine
% computes the tensors \tilde{\hat C}, given by
% \tilde{\hat C}_k = \hat M_L (\hat M_C)^{-1} \hat C_k,
% where \hat M_L and \hat M_C are lumped and consistent mass matrices on
% the reference simplex and \hat C_k = 
% \int_{\hat K} \hat\phi_i \frac{\partial\hat\phi_j}{\partial \hat x_k} dx.
%
% See also H. Hajduk (2021) Monolithic convex limiting in 
% discontinuous Galerkin discretizations of hyperbolic conservation laws 
% Comput. Math. Appl. 87: 120--138 DOI: 10.1016/j.camwa.2021.02.012
%
function CTilde = computePrecMatSimplex(dim, p)
N = nchoosek(p+dim,p);
CTilde = cell(1,dim);
multiindices = getMultiindices(dim,p);
volRefSimplex = 1 / factorial(dim);
M_L = volRefSimplex / N;
gradBary = [-ones(1,dim); eye(dim)];
for m = 1:dim
  CTilde{m} = sparse(N,N);
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
        CTilde{m}(i,j) = CTilde{m}(i,j) + gradBary(k,m)  * (alpha(l) - (l==k) + 1);
      end
    end
  end
  CTilde{m} = M_L * CTilde{m};
end
end
