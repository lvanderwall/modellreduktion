% function [] = MyFun()
%
%   inputs:
%       
%
%   outputs:
%       
% 
%   MyFun(..) ...


function [] = MyFun()
    fprintf('-------- MyFun() --------\n');
    h = waitbar();                  % create waitbar message
    
    % log
    delete(h);                      % delete waitbar
    fprintf('i = %6d,\t e_abs = %.3e\n', i, err(1, i));
    fprintf('-------------------------\n\n');    
end