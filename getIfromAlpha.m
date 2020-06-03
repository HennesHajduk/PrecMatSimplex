% Returns the local index for a given multiindex on the reference triangle
% \hat K = conv\{(0,0),(1,0),(0,1)\} with p-th order polynomials.
% The degrees of freedom are ordered in lexicographical manner.
%
function i = getIfromAlpha(p,alpha)
i = (p+1)*alpha(3) - alpha(3)*(alpha(3)-1)/2 + alpha(2);
i = i+1; % MATLAB counts indices starting at 1.
end