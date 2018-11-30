% function [K] = krylov(A, b, k)
%
%   inputs:
%       A in Mat(n x n, R)
%       b in Mat(n x 1, R)
%       k in N
%
%   outputs:
%       K in Mat(n x k, R)
% 
%   krylov(..) returns the Krylov matrix
%       K(A, b, k) = [b, ..., A^{k-1}b]


function [K] = krylov(A, b, k)
    n = size(A, 2);
    K = zeros(n, k);                % preallocate to improve speed
    
    K(:, 1) = b;
    for j = (2 : 1 : k)             % j > 1: first column already defined
        b = A * b;
        K(:, j) = b;
    end
end