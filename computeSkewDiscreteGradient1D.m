function C = computeSkewDiscreteGradient1D(p)
N = p+1;
C = 0.5 * spdiags([-ones(N,1), ones(N,1)], [-1 1], N, N);
C(1,1) = -0.5;
C(end,end) = 0.5;
end
