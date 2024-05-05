function i = getNodeFromMultiindex(dim, p, alpha)
assert(round(p)==p && p>0, 'Invalid polynomial degree');
assert( isvector(alpha) && length(alpha(:)) == (dim+1) && isequal(round(alpha), alpha) ...
        && min(alpha) >= 0 && sum(alpha)==p, 'Invalid multiindex' );
if dim == 1
  i = alpha(2) + 1;
elseif dim == 2
  i = alpha(2) + 1;
  for j = 1:alpha(3)
    i = i + (p+2-j);
  end
elseif dim == 3
  i = alpha(2) + 1;
  for j = 1:alpha(3)
    i = i + (p+2-j-alpha(4));
  end
  for k = 1:alpha(4)
    i = i + nchoosek(p+3-k,p+1-k);
  end
else
  error('Higher dimensional simplices not supported.');
end
end
