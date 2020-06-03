% Fills an array mapping local indices to multiindeces on the reference 
% triangle \hat K = conv\{(0,0),(1,0),(0,1)\} with p-th order polynomials.
% The degrees of freedom are ordered in lexicographical manner.
%
function i2alpha = mapI2Alpha(p)
dim = 2;                      % spacial dimension
N = nchoosek(p+dim,p);        % number of local degrees of freedom
i2alpha = zeros(N,dim+1);     % initialization
a = [p 0 0];                  % first multiindex
ctr = 1;                      % counter variable
for j = 0:p                   % step through dofs in y-direction
  for i = 0:p-j               % step through dofs in x-direction
    i2alpha(ctr,:) = a;       % fill map entry
    ctr = ctr+1;              % update counter
    a(1:2) = a(1:2) + [-1 1]; % go to next multiindex in x-direction
  end
  a(1) = p-1-j;               % go to next multiindex in y-direction
  a(2) = 0;                   % go to first DOF in x-direction
  a(3) = a(3)+1;              % go to next dof in y-direction
end
end