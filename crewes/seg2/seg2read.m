function varargout=seg2read(sg2file, varargin)
% SEG2READ   Reads SEG-2 format geophysical data (typically Seismic or GPR)
%
% Read all wiggles at once
%   alldata = seg2read(filename);
%
% Read every third trace
%   alldata = seg2read(filename,'traces', [1 3 inf]);  % min, increment, max
%
% Various headers are available if wanted.
%   'data'      2-D matrix (trace,sample) of data [default]
%   'filehdr'   file header (binary) as a structure
%   'filetxt'   file header (text) as a structure
%   'trchdr'    trace headers (binary) as a vector of structures
%   'trctxt'    trace headers (text) as a vector of structures
%
%    [filehdr,filetxt,trchdr,trctxt,data] = seg2read(filename, 'want','filehdr,filetxt,trchdr,trctxt,data');
%
%
% Author: Henry Bland, CREWES, University of Calgary
% Modified: Kevin Hall, CREWES, March 24, 2023
%
% NOTE: This SOFTWARE may be used by any individual or corporation for any purpose
% with the exception of re-selling or re-distributing the SOFTWARE.
% By using this software, you are agreeing to the terms detailed in this software's
% Matlab source file.

% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its 'AUTHOR' (identified above) and the CREWES Project.  The CREWES
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) This SOFTWARE may be used by any individual or corporation for any purpose
%    with the exception of re-selling or re-distributing the SOFTWARE.
%
% 2) The AUTHOR and CREWES must be acknowledged in any resulting publications or
%    presentations
%
% 3) This SOFTWARE is provided "as is" with no warranty of any kind
%    either expressed or implied. CREWES makes no warranties or representation
%    as to its accuracy, completeness, or fitness for any purpose. CREWES
%    is under no obligation to provide support of any kind for this SOFTWARE.
%
% 4) CREWES periodically adds, changes, improves or updates this SOFTWARE without
%    notice. New versions will be made available at www.crewes.org .
%
% 5) Use this SOFTWARE at your own risk.
%
% END TERMS OF USE LICENSE

property_strings = ...
    {'traces','want','verbose'};

if nargin == 0 || isempty(sg2file)
        [name,path]=uigetfile('*.SG2');
        sg2file = fullfile(path,name);
end

if ~ischar(sg2file)
	error('First argument must be a file name');
end

want = FindValue('want',property_strings,varargin{:});
if isempty(want)
    want = 'data';
end



fileinfo = dir(sg2file); % Fileinfo is a structured array of file information.

if isempty(fileinfo)
    error(['The file ', sg2file, ' does not exist.'])
end

fileinfo.bytes; % Pull the size of the file out of fileinfo.
byteOrder = 'le';   % Little-endian byte order

fid = fopen(sg2file, 'r', ['ieee-' byteOrder]);% Open the segy file for reading.

if fid == -1
    error('Unable to open file.')
end

filehdr = readFileDescriptor(fid);
if filehdr.file_descriptor_id ~= hex2dec('3a55')
    error(['Invalid file descriptor: corrupted file? ' sg2file]);
end

ntraces = filehdr.number_of_traces_in_file;

traceDescPtr = readTraceDescriptorPointers(fid, ntraces);

filetxt = readFileParams(fid, filehdr.size_of_trace_pointer_subblock);

% If we want to get traces in some funny order, or skip some,
% we just use a different "wantedTraces" (e.g. 1:3:128)
traces = FindValue('traces',property_strings,varargin{:});
if isempty(traces) 
    wantedTraces = 1:ntraces;
else
    if length(traces) < 2
        traces(2) = 1;
    end
    if length(traces) < 3
        traces(3) = inf;
    end
    wantedTraces = traces(1):traces(2):min(ntraces,traces(3));
end


NoutputTraces = length(wantedTraces);
th = readTraceHeader(fid, traceDescPtr(1)); %read first trace header;
tt = readParamStrings(fid);
% Pre-allocate memory, assuming all traces in file being read are the same
Nsamples = th.number_of_samples;
data = zeros(Nsamples,NoutputTraces);
trchdr = repmat(th,1,NoutputTraces);
trctxt = repmat(tt,1,NoutputTraces);

for seqno=1:NoutputTraces
    traceInFile = wantedTraces(seqno);   % often the same thing as seqno
    th = readTraceHeader(fid, traceDescPtr(traceInFile));
    trchdr(seqno) = th;
    if trchdr(seqno).trace_descriptor_id ~= hex2dec('4422')
       error('crewes:seg2read',['Invalid trace descriptor: corrupted file? ' sg2file]);
    end
    if th.number_of_samples ~= Nsamples
        error('crewes:seg2read','seg2read does not handle variable length traces in a file');
    end
    trctxt(seqno) = readParamStrings(fid);

    advanceToDataBlock(fid, traceDescPtr(traceInFile) +th.size_of_this_block);

    switch trchdr(seqno).data_format_code
        case 1
            data(:,seqno) = fread(fid,trchdr(seqno).number_of_samples,'int16=>double');
        case 2
            data(:,seqno) = fread(fid,trchdr(seqno).number_of_samples,'int32=>double');
        case 3
            error('seg2read does not support 20 bit floating point data');
        case 4
            data(:,seqno) = fread(fid,trchdr(seqno).number_of_samples,'float=>double');
        case 5
            data(:,seqno) = fread(fid,trchdr(seqno).number_of_samples,'double=>double');
        otherwise
            error(['data format ' trchdr(seqno).data_format_code ' is not supported']);
    end
end
fclose(fid);

wanted = split(want,',');
varargout = cell(1,length(wanted)); %pre-allocate cell array

for arg = 1:length(wanted)
    switch wanted{arg}
        case 'data'
            varargout{arg} = data;
        case 'trchdr'
            varargout{arg} = trchdr;
        case 'trctxt'
            varargout{arg} = trctxt;
        case 'filehdr'
            varargout{arg} = filehdr;
        case 'filetxt'
             varargout{arg} = filetxt;
        otherwise
             error('invalid value for "want" parameter. Must be one of: data,trchdr,trctxt,filehdr,filetxt');
    end
end

function advanceToDataBlock(fid, ptr)
    fseek(fid, ptr, 'bof');
    %ftell(fid) %DEBUG
end

function trchdr = readTraceHeader(fid, ptr)
    fseek(fid, ptr, 'bof');
    trchdr.trace_descriptor_id = fread(fid, 1, 'uint16')';
    trchdr.size_of_this_block = fread(fid, 1, 'uint16');
    trchdr.size_of_following_data_block = fread(fid, 1 ,'uint32');
    trchdr.number_of_samples = fread(fid, 1 ,'uint32');
    trchdr.data_format_code = fread(fid, 1, 'uint8');
    trchdr.reserved = fread(fid, 19, 'uint8');
end


function fileParams = readFileParams(fid, size_of_trace_pointer_subblock)
    size_of_file_header = 32;
    fseek(fid, size_of_trace_pointer_subblock + size_of_file_header, 'bof');
    fileParams = readParamStrings(fid);
end

function traceDescPtr = readTraceDescriptorPointers(fid, ntraces)
    traceDescPtr = fread(fid, ntraces, 'uint32');
end

function params = readParamStrings(fid)
    len = fread(fid, 1, 'int16');
    params = struct();
    while len > 0
        s = fread(fid, len - 2, 'char=>char')';
        [name, value] = processParamString(s);
        params.(lower(name)) = value;
        len = fread(fid, 1, 'int16');
    end
end

function [name, value] = processParamString(s)
    [name, value] = strtok(s,' ');
    value(1) = [];
end

function filehdr = readFileDescriptor(fid)
    filehdr.file_descriptor_id = fread(fid, 1, 'uint16')';
    filehdr.revision_number = fread(fid, 1, 'uint16');
    filehdr.size_of_trace_pointer_subblock = fread(fid, 1 ,'uint16');
    filehdr.number_of_traces_in_file = fread(fid, 1 ,'uint16');
    filehdr.size_of_string_terminator = fread(fid, 1, 'uint8');
    filehdr.first_string_terminator = fread(fid, 1, 'uint8');
    filehdr.second_string_terminator = fread(fid, 1, 'uint8');
    filehdr.size_of_line_terminator = fread(fid, 1, 'uint8');
    filehdr.first_line_terminator = fread(fid, 1, 'uint8');
    filehdr.second_line_terminator = fread(fid, 1, 'uint8');
    filehdr.reserved = fread(fid, 18, 'uint8=>uint8');
end

% -----------------------
function value = FindValue( property_name, property_strings, varargin )

    value = [];
    for i = 1:((nargin-2)/2)
        current_name = varargin{2*i-1};
        if ischar(current_name)
            numprops = 1:length(property_strings);
            imatch = numprops(strcmpi(current_name,property_strings));
            nmatch = length(imatch);
            if nmatch > 1
                error(['Ambiguous property name ' current_name '.']);
            end
            if nmatch == 1
                canonical_name = property_strings{imatch};
                if strcmp(canonical_name, property_name)
                    if isempty(value)
                        if isempty(varargin{2*i})
                            error(['Empty value for ' property_name '.']);
                        end
                        value = varargin{2*i};
                    else
                        error(['Property ' property_name ' is specified more than once.']);
                    end
                end
            end
        end
    end
end

end