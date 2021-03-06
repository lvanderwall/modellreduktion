% function [] = my_arnoldi(A, B, C, r, varargin)
%
%   inputs:
%       A in Mat(n x n, R),             system matrix
%       B in Mat(n x m, R),             input matrix
%       C in Mat(p x n, R),             output matrix
%       r in {1; n},                    reduced LTI order
%
%   optional inputs:
%       D = 0 in Mat(p x m, R),         feedthrough matrix
%
%   parameters:
%       in = 1 in {1; m},               input index
%       out = 1 in {1; p},              output index
%       w = logspace(-1, 6, 100)
%           in Mat(1 x len, R)          frequency discretization
%
%
%   my_arnoldi(..) outputs bode- and pole-zero plots of a transfer function
%   of an LTI with input u_in(t), output y_out(t) and state space model
%
%           x^{(1)}(t) = A * x(t) + B * u(t)
%                 y(t) = C * x(t) + D * u(t)
%
%   and its reduced form using the first r arnoldi vectors of Krylov space 
%       K_n(A, b) = Im(b, ..., A^{n-1}b)


function [] = my_arnoldi(A, B, C, r, varargin)
    % ---------------------------------------------------------------------
    % input parsing
    % ---------------------------------------------------------------------
%     n = size(A, 2);                         % dim of state space
    m = size(B, 2);                         % dim of input vector
    p = size(C, 1);                         % dim of output vector

    ip = inputParser;
    ip.FunctionName = 'MY_MODAL';
    
    addRequired(ip, 'A',                            @isnumeric);
    addRequired(ip, 'B',                            @isnumeric);
    addRequired(ip, 'C',                            @isnumeric);
    addRequired(ip, 'r',                            @isnumeric);    % reduced order
    addOptional(ip, 'D',     zeros(p, m),           @isnumeric);
    addParameter(ip, 'in',   1,                     @isnumeric);
    addParameter(ip, 'out',  1,                     @isnumeric);
    addParameter(ip, 'w',    logspace(-1, 6, 100),  @isnumeric);
        
    parse(ip, A, B, C, r, varargin{:});
    
    
    % ---------------------------------------------------------------------
    % reduce LTI using first r arnoldi vectors
    % ---------------------------------------------------------------------
    Q = arnoldi(ip.Results.A, ip.Results.B(:, ip.Results.in), r);

    % calculate reduced LTI
    A_hat = Q' * (ip.Results.A * Q);
    B_hat = Q' * ip.Results.B;
    C_hat = ip.Results.C * Q;
%     D_hat = D;                              % no need to reassign D_hat

    
    % ---------------------------------------------------------------------
    % bode plots
    % ---------------------------------------------------------------------
    n_sub = 3;                              % maximum number of subfigures
    for i = (1 : 1 : n_sub)                 % clear axes
        subplot(n_sub, 1, i);
        cla;
    end
    
    my_bode(...                             % plot unreduced G(iw)
        ip.Results.A,...
        ip.Results.B,...
        ip.Results.C,...
        ip.Results.D,...
        'in',       ip.Results.in,...
        'out',      ip.Results.in,...
        'sub',      n_sub,...
        'w',        ip.Results.w...
    );

    my_bode(...                             % plot reduced G_hat(iw)
        A_hat,...
        B_hat,...
        C_hat,...
        ip.Results.D,...                    % = D_hat
        'in',       ip.Results.in,...
        'out',      ip.Results.in,...
        'sub',      n_sub,...
        'w',        ip.Results.w...
    );

    % ---------------------------------------------------------------------
    % pole-zero plot
    % ---------------------------------------------------------------------
    S_un = eigs(A, size(A, 2));             % eigenvalues of unreduced LTI
    S = eig(A_hat);                         % eigenvalues of reduced LTI
    
    subplot(n_sub, 1, 3);                   % my_bode(..) uses subplot(n_sub, 1, 1-2) by default
    hold off;
    plot(...
        real(diag(S_un)), imag(diag(S_un)), 'ok',...
        real(diag(S)),    imag(diag(S)),    '*r'...
    );
    xlabel('Re{}');
    ylabel('Im{}');
    title(sprintf('\\lambda_k of G_{%i,%i}(i\\omega) and G_{%i,%i,hat}(i\\omega)',...
        ip.Results.out, ip.Results.in,...
        ip.Results.out, ip.Results.in...
    ));
end