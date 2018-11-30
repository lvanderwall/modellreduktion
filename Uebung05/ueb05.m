% -------------------------- Modellreduktion A05 --------------------------
% Name: Lennert van der Wall
% -------------------------------------------------------------------------
%
%
% table of contents
%   {A1}
%   {A2}

dbstop in ueb05.m at 51;



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



% -------------------------------------------------------------------------
% {A2}
% -------------------------------------------------------------------------
% sites:
%   http://slicot.org/20-site/126-benchmark-examples-for-model-reduction
%   
% 
% files: 
%     http://slicot.org/objects/software/shared/bench-data/CDplayer.zip
clear variables;
clc;

load CDplayer.mat;

w = logspace(-1, 6, 300);           % frequency discretization in s^{-1}
r = (5 : 1 : 118);                  % reduction orders to try

h = waitbar(0, 'reduction order k =   5...');
figure;
for k = r                           % compare reduced to unreduced model
    my_arnoldi(...
        A, B, C, k,...
        'w', w...
    );
    pause(0.1);
    waitbar(k / r(end), h, sprintf('reduction order k = %3i...', k));
end
delete(h);


% Compared to HA3_A1), the eigenvalues of the reduced system are NOT a
% subset of \sigma(A). That's not surprising, as we use orthonormal
% projection to a krylov subspace K_r(A, b) that is generally NOT 
% A-invariant, so it generally is not a direct sum of eigen- and hauptspaces
% of A.
%
% For high frequency ranges, the reduced system approximates the amplitude 
% function well with a decreasing lower bound when rank r increases.
%
% Complex conjugate eigenpairs should NOT be approximated seperately:
% Phase- and amplitude approximation seem to get worse if one real 
% eigenvalue approximates a complex pair!