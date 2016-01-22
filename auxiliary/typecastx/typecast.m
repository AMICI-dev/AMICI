%  Y = typecast(X, type) converts a numeric value in X to the data type specified by type.
%*************************************************************************************
% 
%  MATLAB (R) is a trademark of The Mathworks (R) Corporation
% 
%  Function:    typecast
%  Filename:    typecast.m
%  Programmer:  James Tursa
%  Version:     3.00
%  Date:        March 17, 2011
%  Copyright:   (c) 2009, 2011 by James Tursa, All Rights Reserved
%
%  This code uses the BSD License:
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.
%
% typecast is a mex function intended to mimic the MATLAB intrinsic typecast function
% for those users with older versions of MATLAB that do not have this intrinsic. Users
% of newer versions of MATLAB may be interested in this C-mex version to take advantage
% of the several extensions offered. This C-mex version of typecast differs from the
% intrinsic typecast in the following important aspects: 
% 
%                      Intrinsic typecast    C-mex typecast
%                      ------------------    --------------
% Type of copy:        Deep Data Copy        Shared Data Copy
% Allows complex:      No                    Yes
% Allows logical:      No                    Yes (cannot convert from complex)
% Allows char:         No                    Yes (cannot convert from complex)
% Allows non-vector:   No                    Yes
%
% Since this C-mex typecast produces a shared data copy of the original, it
% is more efficient than the MATLAB intrinsic typecast, which may be important
% if you are working with large variables. For non-vector inputs, the first
% non-singleton dimension must be compatible for the conversion.
%
% Building:
%
% typecast requires that a mex routine be built (one time only). This
% process is typically self-building the first time you call the function
% as long as you have the files typecast.m and typecast.c in the same
% directory somewhere on the MATLAB path. If you need to manually build
% the mex function, here are the commands:
%
% >> mex -setup
%   (then follow instructions to select a C or C++ compiler of your choice)
% >> mex typecast.c
%
% The usage is as follows (from The Mathworks website documentation):
%
% Syntax
%
% Y = typecast(X, type)
%
% Description
%
% Y = typecast(X, type) converts a value in X to the data type specified by type.
% Input X must be a full numeric or char or logical variable. The type input is a string
% set to one of the following: 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
% 'uint64', 'int64', 'single', 'double', 'char' or 'logical'. typecast is different
% from the MATLAB cast function in that it does not alter the input data. typecast
% always returns the same number of bytes in the output Y as were in the input X.
% For example, casting the 16-bit integer 1000 to uint8 with typecast returns the full
% 16 bits in two 8-bit segments (3 and 232) thus keeping its original value
% (3*256 + 232 = 1000). The cast function, on the other hand, truncates the input value
% to 255.
%
% The output of typecast can be formatted differently depending on what system you use it on.
% Some computer systems store data starting with its most significant byte (an ordering
% called big-endian), while others start with the least significant byte (called little-endian). 
%
% typecast issues an error if X contains fewer values than are needed to make an output value. 
% 
%*************************************************************************************
 
function varargout = typecast(varargin)
fname = 'typecast';
disp(' ');
disp(['Detected that the mex routine for ' fname ' is not yet built.']);
disp('Attempting to do so now ...');
disp(' ');
try
    mname = mfilename('fullpath');
catch
    mname = fname;
end
cname = [mname '.c'];
if( isempty(dir(cname)) )
    disp(['Cannot find the file ' fname '.c in the same directory as the']);
    disp(['file ' fname '.m. Please ensure that they are in the same']);
    disp('directory and try again. The following file was not found:');
    disp(' ');
    disp(cname);
    disp(' ');
    error(['Unable to compile ' fname '.c']);
else
    disp(['Found file ' fname '.c in ' cname]);
    disp(' ');
    disp('Now attempting to compile ...');
    disp('(If prompted, please press the Enter key and then select any C/C++');
    disp('compiler that is available, such as lcc.)');
    disp(' ');
    disp(['mex(''' cname ''',''-output'',''',mname,''')']);
    disp(' ');
    try
        mex(cname,'-output',mname);
        disp([ fname ' mex build completed ... you may now use ' fname '.']);
        disp(' ');
    catch
        disp(' ');
        error(['Unable to compile ' fname ' ... Contact author.']);
    end
    [varargout{1:nargout}] = typecast(varargin{:});
end
end
