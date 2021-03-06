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
N = nchoosek(p+dim,p);    % amount of local DOF
CTilde = zeros(N,N,dim);  % initialization
i2Alpha = mapI2Alpha(p);  % array mapping local indices to multiindices
volRefSimplex = 1 / factorial(dim); % |\hat K|
M_L = volRefSimplex / N;  % diagonal entries of lumped mass matrix \hat M_L
gradBary = [-ones(1,dim); eye(dim)]; % gradients of barycentric coordinates
for m = 1:dim             % calc. \tilde{ \hat C}_m
  for j = 1:N             % calc. (\tilde{ \hat C}_m)_{:,j}
    alpha = i2Alpha(j,:); % multiindex for jth local index
    for k = 1:dim+1
      for l = 1:dim+1
        beta = alpha;
        beta(k) = beta(k) - 1;
        beta(l) = beta(l) + 1;
        if (beta(l) > p || beta(k) < 0)
          continue;
        end
        i = getIfromAlpha(p,beta); % local index for beta
        CTilde(i,j,m) = CTilde(i,j,m) ...   % cf. % (B.1)
            + gradBary(k,m)  * (alpha(l) - (l==k) + 1);
      end
    end
  end
end
CTilde = M_L * CTilde; % multiply with lumped mass matrix
end
