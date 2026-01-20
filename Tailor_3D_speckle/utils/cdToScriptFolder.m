function cdToScriptFolder()
% cdToScriptFolder  Change current folder to the folder where the calling script resides
%
% Usage:
%   Call this at the top of any script or function to switch MATLAB's
%   working directory to the directory containing that file.

    % Get the call stack with full file names
    stk = dbstack('-completenames');
    if numel(stk) < 2
        % If called from the command line or anonymous function, nothing to do
        return
    end

    % The caller is the second entry in the stack
    callerFile = stk(2).file;            % full path to the .m file that called this
    callerFolder = fileparts(callerFile);

    % Change to that folder
    cd(callerFolder);
end
