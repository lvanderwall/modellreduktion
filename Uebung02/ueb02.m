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


% -------------------------------------------------------------------------
% {A1b}
% -------------------------------------------------------------------------
% plot amplitudes, cannot use my_bode to get everything in one plot:
w = logspace(-1, 2, 1e3);                   % w in s^{-1}

% plot reduced models
figure;                                     % new figure
for k = (1 : 1 : 3)
    tf = @(s) (...                          % scalar G_k,hat(iw)
        C_hat{k} * (...
            (s .* eye(size(A_hat{k}, 2)) - A_hat{k}) \ B_hat{k}...
        )...
    );

    G = arrayfun(tf, 1i * w);               % calculate complex tf values
    semilogx(w, 20 * log10(abs(G)));
    hold on;
end

% plot unreduced model
tf = @(s) (...                              % scalar G(iw)
    c * (...
        (s .* eye(size(A, 2)) - A) \ b...
    )...
);

G = arrayfun(tf, 1i * w);                   % calculate complex tf values
semilogx(w, 20 * log10(abs(G)), '-k');      % bode plot of unreduced LTI

% annotations
title("bode plot for G_{1, 1}(i\omega)");
xlabel("\omega in s^{-1}");
ylabel("|G_{1, 1}(i\omega)| in dB20");
