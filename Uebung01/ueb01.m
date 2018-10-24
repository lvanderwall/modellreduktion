% -------------------------- Modellreduktion A ----------------------------
% Name: Lennert van der Wall
% -------------------------------------------------------------------------
%
%
% table of contents
%   {A1}
%   {A2}

dbstop in ueb01.m at 41; dbstop in ueb01.m at 73;



% -------------------------------------------------------------------------
% {A1}
% -------------------------------------------------------------------------
% sites:
%   http://slicot.org/20-site/126-benchmark-examples-for-model-reduction
%   
% 
% files: 
%     http://slicot.org/objects/software/shared/bench-data/CDplayer.zip
%     https://sparse.tamu.edu/mat/Oberwolfach/rail_1357.mat
clear variables;
clc;

% -------------------------------------------------------------------------
% a)
% -------------------------------------------------------------------------
fprintf('---------------- A1a -----------------\n');
load CDplayer.mat;

w = logspace(-1, 6, 300);           % frequency discretization in s^{-1}
my_bode(A, B, C, 'w', w);
my_sigma(A, B, C, 'w', w);

% -------------------------------------------------------------------------
% b)
% -------------------------------------------------------------------------
fprintf('---------------- A1b -----------------\n');
clear variables;
clc;

load rail_1357.mat;
m = size(Problem.aux.B, 2);         % size of input vector
p = size(Problem.aux.C, 1);         % size of output vector
D = zeros(p, m);                    % D matrix missing

w = logspace(-6, 2, 300);           % frequency discretization in s^{-1}
my_bode(...
    Problem.A,...
    Problem.aux.B,...
    Problem.aux.C,...
    D,...                 
    Problem.aux.E,...
    'w', w...
 );

my_sigma(...
    Problem.A,...
    Problem.aux.B,...
    Problem.aux.C,...
    D,...
    Problem.aux.E,...
    'w', w...
 );

 
 
% -------------------------------------------------------------------------
% {A2}
% -------------------------------------------------------------------------
clear variables;
clc;

for k = (60 : 20 : 110)
	my_img_SVD('einstein.jpg', k);
end

% k = 60 -> 110:
%	notice artifacts only in background and hair
%
% k > 110:
%	artifacts in background nearly invisible, only hair remains blurred