% -------------------------- Modellreduktion A05 --------------------------
% Name: Lennert van der Wall
% -------------------------------------------------------------------------
%
%
% table of contents
%   {A1}
%   {A2}

% dbstop in ueb_def.m at 1;



% -------------------------------------------------------------------------
% {A1}
% -------------------------------------------------------------------------
clear variables;
clc;


% -------------------------------------------------------------------------
% test system

A = [...
    1, 2, 3;...
    4, 5, 6;...
    7, 8, 9 ...
];

k = 3;

b = ones(k, 1);

% -------------------------------------------------------------------------
% test arnoldi(..) and krylov(..)
Q = arnoldi(A, b, k);               % should return dim(Q) = 1 < k
if(size(Q, 2) == 1)
    disp('arnoldi(..) dimension < k tested successfully!');
end


b = [1, 1, 2]';
Q = arnoldi(A, b, k);                % should return dim(Q) = 3 = k
if(size(Q, 2) == 3)
    disp('arnoldi(..) dimension = k tested successfully!');
end

K = krylov(A, Q(:, 1), k);
H = Q' * A * Q;
R = krylov(H, eye(3, 1), k);
if(norm(K - Q * R, 'fro') < 1e-12)
    disp('K(q1, A, k) = Q * R with R = K(Q^TAQ, e1, k) tested successfully!');
end

