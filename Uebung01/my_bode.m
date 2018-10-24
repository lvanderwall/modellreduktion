% function [] = my_bode(A, B, C, varargin)
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
%   parameters:
%       in = 1 in {1; m},               input index
%       out = 1 in {1; p},              output index
%       w = logspace(-1, 6, 100)
%           in Mat(1 x len, R)          frequency discretization
%
%
%   my_bode(..) plots a bode plot of a transfer function of a generalized
%   state space with input u_in(t), output y_out(t) and state space model
%
%       E * x^{(1)}(t) = A * x(t) + B * u(t)
%                 y(t) = C * x(t) + D * u(t)
%
%   It plots |G(iw)|_{dB20} and \phi(G(iw) over log10(w)


function my_bode(A, B, C, varargin)
    % ---------------------------------------------------------------------
    % input parsing
    % ---------------------------------------------------------------------
    n = size(A, 2);                         % dim of state space
    m = size(B, 2);                         % dim of input vector
    p = size(C, 1);                         % dim of output vector

    ip = inputParser;
    ip.FunctionName = 'MY_BODE';
    
    addRequired(ip, 'A',                            @isnumeric);
    addRequired(ip, 'B',                            @isnumeric);
    addRequired(ip, 'C',                            @isnumeric);
    addOptional(ip, 'D',     zeros(p, m),           @isnumeric);
    addOptional(ip, 'E',     speye(n),              @isnumeric);    % use speye to optimize memory usage
    addParameter(ip, 'in',   1,                     @isnumeric);
    addParameter(ip, 'out',  1,                     @isnumeric);
    addParameter(ip, 'w',    logspace(-1, 6, 100),  @isnumeric);
        
    parse(ip, A, B, C, varargin{:});
    
    
    % ---------------------------------------------------------------------
    % plotting
    % ---------------------------------------------------------------------
    tf = @(s) (...                          % scalar G_{out, in}(iw)
        ip.Results.C(ip.Results.out, :) * (...
            (s .* ip.Results.E - ip.Results.A) \ ip.Results.B(:, ip.Results.in)...
        ) + ip.Results.D(ip.Results.out, ip.Results.in)...
    );

    G = arrayfun(tf, 1i * ip.Results.w);    % calculate complex tf values
    
    figure;                                 % new figure
    subplot(2, 1, 1);
    semilogx(ip.Results.w, 20 * log10(abs(G)));
    title(sprintf("bode plot for G_{%i, %i}(i\\omega)", ip.Results.out, ip.Results.in));
    xlabel("\omega in s^{-1}");
    ylabel(sprintf("|G_{%i, %i}(i\\omega)| in dB20", ip.Results.out, ip.Results.in));

    subplot(2, 1, 2);
    semilogx(ip.Results.w, angle(G) / pi);
    title(sprintf("bode plot for G_{%i, %i}(i\\omega)", ip.Results.out, ip.Results.in));
    xlabel("\omega in s^{-1}");
    ylabel(sprintf("\\angle G_{%i, %i}(i\\omega) in \\pi", ip.Results.out, ip.Results.in));  
end