% Fills an array mapping local indices to multiindeces on the reference 
% triangle \hat K = conv\{(0,0),(1,0),(0,1)\} with p-th order polynomials.
% The degrees of freedom are ordered in lexicographical manner.
%
function i2alpha = mapI2Alpha(p)
dim = 2;                      % spatial dimension
N = nchoosek(p+dim,p);        % amount of local DOF
i2alpha = zeros(N,dim+1);     % initialization
a = [p 0 0];                  % first DOF
ctr = 1;                      % counter variable
for j = 0:p                   % step through DOF in y-dir.
  for i = 0:p-j               % step through DOF in x-dir.
    i2alpha(ctr,:) = a;       % fill map entry
    ctr = ctr+1;              % update counter
    a(1:2) = a(1:2) + [-1 1]; % go to next DOF in x-dir.
  end
  a(1) = p-1-j;               % go to next DOF in y-dir.
  a(2) = 0;                   % go to first DOF in x-dir.
  a(3) = a(3)+1;              % go to next DOF in y-dir.
end
end
