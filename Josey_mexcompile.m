function okay = Josey_mexcompile(funcName, varargin)
%VGG_MEXCOMPILE  Mex compile helper function
%
%   okay = vgg_mexcompile(funcName)
%   okay = vgg_mexcompile(..., sourceList)
%   okay = vgg_mexcompile(..., lastCompiled)
%
% Compile mexed function, given an optional list of source files. Can
% optionally check if source files have been modified since the last
% compilation, and only compile if they have.
%
%IN:
%   funcName - string containg the name of the function to compile
%   sourceList - cell array of source files to be compiled. Default:
%                {[funcName '.cxx']}.
%   lastCompiled - datenum of the current mex file. Default: 0 (i.e. force
%                  compilation).
%
%OUT:
%   okay - 1: function compiled; 0: up-to-date, no need to compile; -1:
%          compilation failed.

% $Id: vgg_mexcompile.m,v 1.3 2009/09/13 20:34:58 ojw Exp $

% Set defaults for optional inputs

% Parse inputs
sourceList = cell(2,1);
sourceList{1} = 'TrainGMM.cpp';
sourceList{2} = 'PredictGMM.cpp';
sourceDir = 'C:\opencv2.3\build\include';
libDir = 'C:\opencv2.3\build\x64\vc10\lib';
libList = '-lopencv_core243 -lopencv_ml243 -lopencv_highgui243';
% Compile if we need to

% Set the compiler flags
flags = ['-O -I' sourceDir ' -L' libDir ' ' libList ];
switch mexext
    case 'mexsol'
        flags = [flags ' CC=gcc CFLAGS=-fPIC'];
    case {'mexglx', 'mexa64'}
        str = '"-O3 -ffast-math -funroll-loops"';
        flags = sprintf('%s CXXOPTIMFLAGS=%s LDCXXOPTIMFLAGS=%s LDOPTIMFLAGS=%s', flags, str, str, str);
    otherwise
end

% Call mex to compile the code
for i=1:size(sourceList,1)
    cmd = sprintf('mex %s%s', flags, sprintf(' "%s"', sourceList{i}));
    disp(cmd);
    try
        eval(cmd);
        okay = 1;
    catch
        okay = -1;
        fprintf('ERROR while compiling %s\n', funcName);
    end

end
% Return to the original directory

return
