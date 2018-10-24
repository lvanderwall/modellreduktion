% function [] = my_sigma(A, B, C, varargin)
%
%   inputs:
%       A in Mat(n x n, R),             system matrix
%       B in Mat(n x m, R),             input matrix
%       C in Mat(p x n, R),             output matrix
%
%   optional inputs:
%       D = 0 in Mat(p x m, R),         feedthrough matrix
%       E = 1 in Mat(n x n, R),         description matrix
%
%
%   my_sigma(..) plots a sigma plot of a transfer function of a generalized
%   state space with input u(t), output y(t) and state space model
%
%       E * x^{(1)}(t) = A * x(t) + B * u(t)
%                 y(t) = C * x(t) + D * u(t)
%
%   It plots log10(sigma_max(w)) and log10(sigma_min(w)) over log10(w)


function my_sigma(A, B, C, varargin)
    % ---------------------------------------------------------------------
    % input parsing
    % ---------------------------------------------------------------------
    n = size(A, 2);                 % dim of state space
    m = size(B, 2);                 % dim of input vector
    p = size(C, 1);                 % dim of output vector

    ip = inputParser;
    ip.FunctionName = 'MY_BODE';
    
    addRequired(ip, 'A',                            @isnumeric);
    addRequired(ip, 'B',                            @isnumeric);
    addRequired(ip, 'C',                            @isnumeric);
    addOptional(ip, 'D',     zeros(p, m),           @isnumeric);
    addOptional(ip, 'E',     speye(n),              @isnumeric);
    addParameter(ip, 'w',    logspace(-1, 6, 100),  @isnumeric);
        
    parse(ip, A, B, C, varargin{:});
    
    
    % ---------------------------------------------------------------------
    % plotting
    % ---------------------------------------------------------------------
    tf = @(s) (...                  % matrix valued G_(iw)
        ip.Results.C * (...
            (s .* ip.Results.E - ip.Results.A) \ ip.Results.B...
        ) + ip.Results.D...
    );
    
    % calculate all sigma(w)
    Sigma = arrayfun(@(s) svd(tf(s)), 1i * ip.Results.w, 'UniformOutput', false);

    figure;                         % new figure
    loglog(...
        ip.Results.w, cellfun(@(A) A(1),    Sigma), '-r',...    % sigma_max(w)
        ip.Results.w, cellfun(@(A) A(end),  Sigma), '-k'...     % sigma_min(w)
    );
    title("sigma plot for G(i\omega)");
    xlabel("\omega in s^{-1}");
    ylabel("\sigma_{min}(\omega), \sigma_{max}(\omega)");
end