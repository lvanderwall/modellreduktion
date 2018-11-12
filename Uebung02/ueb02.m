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



% -------------------------------------------------------------------------
% {A2}
% -------------------------------------------------------------------------
clear variables;
open 'ww_36_pmec_36.mat';           % SISO LTI without E of order n = 66

% rename matrizes, do not use ans!
A = ans.A;                          % in Mat(66 x 66, R)
b = ans.b;                          % in Mat(66 x  1, R)
c = ans.c;                          % in Mat( 1 x 66, R)
d = ans.d;                          % in R
clc;

% -------------------------------------------------------------------------
% {A2a}
% -------------------------------------------------------------------------
n = size(A, 2);
g_amp = @(s) (...                   % |G(s)|_dB20
    20.0 * log10(abs(...
        c * (...
            (s * eye(n) - A) \ b...
        ) + d...
    ))...
);

points = 200;                       % number of points per axis
re_data = linspace(-2, 0, points);  % real axis
im_data = linspace(0, 20, points);  % imaginary axis

% amplitudes in db20, preallocate to improve speed
am_data = zeros(size(re_data, 2), size(im_data, 2));
h = waitbar(0, 'Calculate amplitudes...');
for i = (1 : 1 : size(re_data, 2))
    for j = (1 : 1 : size(im_data, 2))
        am_data(i, j) = g_amp(re_data(i) + 1.0i*im_data(j));
    end
    waitbar(i / size(re_data, 2), h);
end
delete(h);

% get poles and index of dominant poles
poles = eig(A);
[~, pole_idx] = sort(real(poles), 'descend');


% im -> x-axis, re -> y-axis, am -> z-axis
figure;                             % new figure
h = meshc(im_data, re_data, am_data);
hold on;
plot3(...                           % |G(iw)|_dB20
    im_data,...
    zeros(size(im_data)),...
    am_data(end, :),...
    '-k',...
    'LineWidth', 2 ...
);

% poles
z_data = min(h(1).Parent.ZAxis.Limits);
plot3(...
    imag(poles),...
    real(poles),...
    z_data*ones(size(poles)),...
    'ok',...
    'MarkerFaceColor', [0, 0, 0]...
);

% dominant poles, there seems to be one instable pole at 8e-5!
n_dom = 10;                         % only 10 dominant poles
plot3(...
    imag(poles(pole_idx(1 : 1 : n_dom))),...
    real(poles(pole_idx(1 : 1 : n_dom))),...
    z_data*ones(1, n_dom),...
    'ok',...
    'MarkerFaceColor', [0, 1, 0]...
);

% annotation
h(1).Parent.YDir = 'reverse';       % reverse y-axis
xlabel('Im(s)');
ylabel('Re(s)');
zlabel('|G(s)|_{dB20}');
xlim([min(im_data), max(im_data)]); % do not expand axes for poles plot
ylim([min(re_data), max(re_data)]);
