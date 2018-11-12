% -------------------------- Modellreduktion A02 --------------------------
% Name: Lennert van der Wall
% -------------------------------------------------------------------------
%
%
% table of contents
%   {A1}
%       {A1a}
%       {A1b}
%   {A2}

dbstop in ueb02.m at 98; dbstop in ueb02.m at 121;


% -------------------------------------------------------------------------
% {A1}
% -------------------------------------------------------------------------
clc;
clear variables;
close all;

% LTI
A = [...
    -2.5,   -.5,    0;...
    -.5,    -2.5,   0;...
    -8.0    4.0,    -5 ...
];
b = [2.0; 0; 5.0];                  % column vector
c = [29.5,  -6.5,   7.0];


% -------------------------------------------------------------------------
% {A1a}
% -------------------------------------------------------------------------
% define projection components and calculate reduced model matrices.
% cells can have any objects as elements:
S = {...                            % left projection components
        [...                        % S_1
            1,  0;...
            1,  0;...
            -2, 1 ...
        ],...        
        [...                        % S_2
            -1, 0;...
            1,  0;...
            4,  1 ...
        ],...
        [...                        % S_3
            -1, 1;...
            1,  1;...
            4,  -2 ...
        ]...
    };

W = {...                      % right projection components
        1 / 2 * [...                % W_1
            1,  6;...
            1,  -2;...
            0,  2 ...
        ],...
        1 / 2 * [...                % W_2
            -1, 6;...
            1,  -2;...
            0,  2 ...
        ],...
        1 / 2 * [...                % W_3
            -1, 1;...
            1,  1;...
            0,  0 ...
        ]...
    };

% define reduced models
for k = (1 : 1 : 3)
    A_hat{k} = W{k}' * A * S{k};
    B_hat{k} = W{k}' * b;
    C_hat{k} = c * S{k};
end