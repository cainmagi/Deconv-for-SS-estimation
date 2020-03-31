function [B] = create_spline_mat(n_y, delta, delta_knot)
% create_spline_mat
%   Use this function to create spline matrix. This function is extracted
%   from the cvxEDA.
%   Arguments:
%       n_y: length of the signal y (required).
%       delta: sampling interval (in seconds) of y. (default 1).
%       delta_knot: time between knots of the tonic spline function
%                   (default 10).
if nargin < 2
    delta = 0.01;
end
if nargin < 3
    delta_knot = 10;
end
delta_knot_s = round(delta_knot / delta);
spl = [1:delta_knot_s delta_knot_s-1:-1:1]'; % order 1
spl = conv(spl, spl, 'full');
spl = spl / max(spl);
% matrix of spline regressors
i = bsxfun(@plus, (0:length(spl)-1)'-floor(length(spl)/2), 1:delta_knot_s:n_y);
nB = size(i, 2);
j = repmat(1:nB, length(spl), 1);
p = repmat(spl(:), 1, nB);
valid = i >= 1 & i <= n_y;
B = full(sparse(i(valid), j(valid), p(valid)));

% trend
C = [ones(n_y,1), (1:n_y)'/n_y];
B = [B, C];
end

