function C = computePrecMat1D(p)
n = p+1;
C = { spdiags([-1:-1:-n; 2*(0:p)-p; p+1:-1:1]', -1:1, n, n) / n };
end
