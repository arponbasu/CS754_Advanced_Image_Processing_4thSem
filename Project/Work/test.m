A = 2*binornd(1,0.5,100,500)-1;
x_orig = binornd(1,0.1,500,1);
y = A*x_orig;
max_iter= 100;
epsilon = 0.01;

% Reweighted L1 Minimization
% min |x|_1
% s.t. Ax=b
%
% INPUT 
%   A: Matrix n x p where p > n
%   b: Vector n x 1
%   max_iter: Max number of iterations to run
% 
% Paper link: 
%   https://web.stanford.edu/~boyd/papers/pdf/rwl1.pdf
[~, p] = size(A);
w = ones(p, 1);
w_old = w;

for i = 1:max_iter
  
  W = diag(w);
  
  cvx_begin quiet
    variable x(p,1)
    minimize(norm(W*x, 1))
    subject to
      A*x == y;
  cvx_end
  val = sum(abs(x)>1e-3);
  disp(val);
  diff = sum(abs(A*x-y));
  disp(diff);
  
  w = 1./(epsilon + abs(x));
  if norm(w-w_old, 2) < 1e-9
    break
  end
  w_old = w;
end

x(abs(x)<1e-3) = 0;
final_val=sum(abs(x_orig)>0);
fprintf('\nx-orig-')
disp(final_val);
fprintf('\ndifference-')
final=sum(abs(x-x_orig)>1e-3);
disp(final);



