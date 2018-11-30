% -------------------------- Modellreduktion A03 --------------------------
% Name: Lennert van der Wall
% -------------------------------------------------------------------------
%
%
% table of contents
%   {A1}



% -------------------------------------------------------------------------
% {A1}
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
    my_modal(...
        A, B, C, k,...
        'w', w,...
        'sigma', 'smallestabs'...
    );
    pause(0.1);
    waitbar(k / r(end), h, sprintf('reduction order k = %3i...', k));
end
delete(h);

% possible strategies
%     'smallestabs':    |G(iw)| seems to be better approximated for small w
%     'largestabs':     arg(Giw) seems to be well approximated for large w
%
%     'largestreal':    those three do not converge well
%     'smallestreal':
%     'bothendsreal':
%
% convergence problems for the last three strategies might arise because
% there are a lot of eigenvalues clustered around a specific real value.
% MATLAB doc of eigs(..) hints that these strategies use a krylov method 
% that seems to struggle when eigenvalues that have a similar real value
%
% I am NOT sure about that reasoning as I do not know enough about the
% methods used!
%
% Anyway, I guess one should NOT seperate a complex conjugate pole pair - 
% that seems to greatly warp the argument of the bode plot!