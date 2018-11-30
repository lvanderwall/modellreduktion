% function [Q] = arnoldi(A, b, k)
%
%   inputs:
%       A in Mat(n x n, R)
%       b in Mat(n x 1, R)
%       k in N
%
%   outputs:
%       Q in Mat(n x l, R) with l <= min(n, k)
% 
%   arnoldi(A, b, k) returns an orthonormal basis of the Krylov subspace
%           K_k(A, b) = im([b, ..., A^{k-1}b])
%   with dim(K_k(A, b)) = l

%#ok<*SPRIX>                        % sparse matrices are preallocated, suppress corresponding warnings

function [Q] = arnoldi(A, b, k)
    n = size(A, 2);                 % size of square matrix A
    eps = 1e-12;                    % ||q||_2 < eps => q ~ 0
    
    H = spalloc(n, n,...            % preallocate to improve speed
        n*(n + 3)/2 - 1);           % H is hessenberg matrix => sparse matrix
        
    Q = zeros(n, k);                % preallocate to improve speed
    
    % arnoldi algorithm
    Q(:, 1) = b / norm(b, 2);
    for j = (1 : 1 : (k - 1))       % K_k(A, b) has only (k-1) + 1 basis vectors!
        z = A * Q(:, j);
        
        for i = (1 : 1 : j)
            H(i, j) = Q(:, i)' * z;
            z = z - H(i, j) * Q(:, i);  % Q is ON - matrix => q_i^T * A * q_j = q_i^T * z
        end
        
        rj = norm(z, 2);            % ||z||_2 < eps => K_{j + 1}(A, b) ~ K_j(A, b)
        if(rj < eps)
            Q = Q(:, j);
            break;
        end
        
        H(j + 1, j) = rj;
        Q(:, j + 1) = z / rj;
    end
end