function [x, mu] = quadprog_gcv(y, A, mu, mu_max)
%quadprog_gcv Quadratic programming with GCV
%   Solve the inverse problem of
%       ||y - Ax||^2_2 + mu || x ||^2_2, where x>=0, Ax <= y
%   with GCV. GCV is used for determining the value of mu.
%       The plain problem could be rewritten as
%       
%   Arguments:
%       y: the target vector
%       A: the dictionary matrix
%       mu_max: maximal value of the mu.
%   Returns:
%       x: the inverted signal.
options = optimoptions('quadprog', 'Display','off');
n_x = size(A,2);
H = A'*A + mu * eye(n_x);
f = -y'*A;
x = quadprog(H, f, A, y, [], [], [zeros([n_x-2,1]); zeros(2,1)], [], [], options);
if nargout > 1
    mu = fminbnd(@(m) gcv_minmizer(m, A, y), 0, mu_max);
end
end

function gcv = gcv_minmizer(mu, A, y)
    Ax = A(:,1:(end-2));
    n_x = size(Ax,2);
    H = Ax / (Ax' * Ax + mu * eye(n_x)) * Ax';
    H_tr = trace(H);
    gcv = mean((H*y).^2, 'all')/H_tr;
end