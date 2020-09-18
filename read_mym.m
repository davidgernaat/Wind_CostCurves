function [data,time,header,error,msg] = read_mym(filename)
% [data,time,header,error,msg] = read_mym(filename)
% purpose:  read a MyM data file (with/without time variable), assuming
%           single variable per file.
%           <header>    contains the leading header lines
%                       where the last line ends with [a1...an](t) = [
%           <data>      contains the data as a [a1...an](t) matrix
%           <time>      vector with time entries (length t)
%           <error>     error indicator (0 = no error)
%           <msg>       error message (empty on succes)
%
% version:  1.0 / 20150608 / MvdB

%% INITIALIZATION & OPEN FILE

% Initialize output variables
error    = 0;
msg      = [];
data     = [];
time     = [];
header   = [];

% Initialize intermediate variables
dataTime = [];  % data variable for single timestep
timeVar  = 1;   % indicates whether variable is time-dependent
                % <timeVar> = 0 means time-independent
n        = 0;   % counter for #data entries
t        = 0;   % counter for #time entries

% Check input arguments
if nargin ~= 1
    error = 1;
    msg   = '1 input variable required';
elseif ~ischar(filename)
    error = 2;
    msg   = 'input variable must be a valid string';
elseif ~exist(filename, 'file') == 2
    error = 3;
    msg   = 'input variable must be an existing file';
end

if error
    return;
end

% Open file
fid = fopen(filename,'rt');
if fid < 0
    error = 4;
    msg   = 'error opening file';
    return;
end

% Read header information
line   = fgetl(fid);            %read line
line   = strtrim(line);         %remove leading/trailing blanks
while strcmp(line(1),'!')
    header = [header;{line}];
    line = fgetl(fid);          %read line
    line = strtrim(line);       %remove leading/trailing blanks
end

%% DETERMINE VARIABLE DIMENSIONS

% Read dimensional information 
% the next line contains "[a1,a2,...,an](t) = [", with ai n dimensions
% the second "[" is optional and can be followed by the first time value

% Determine whether file contains time-dependent variable
indtime=strfind(line,'(t)');
if isempty(indtime),
    timeVar = 0;
end

ind1=strfind(line,'[');
ind2=strfind(line,']');

% Extract possible first year value after '=' or 2nd '['
if line(end) ~= '[' && line(end) ~= '='
    firstYear = str2double(regexp(line(ind1(end):end), '[0-9]*', 'match'));
    dataTime  = firstYear;
    n         = 1;
end

% Extract dimensional information, place in variable <dims>
if isempty(ind2)
    msg = 'No dimensional information';
    dims = [1,1];
else
    line = line(ind1(1):ind2(1));
    dims = str2double(regexp(line, '[0-9]*', 'match'));
end

%% READ DATA

% Create data array with size(dims)
nPerTime = prod(dims);
data     = zeros(1,nPerTime);

% Read in data
while ~feof(fid)
    % read in data until the required size <nPerTime> is reached    
    while n < nPerTime + timeVar && ~feof(fid)
        % read line
        line = fgetl(fid);
        % transform data to double and update counter <n>, using
        % regular expressions, for either normal or scientific notation.
        dataLine = str2double(regexp(line, '[0-9.eE+\-]*', 'match'));
        nLine    = numel(dataLine);
        % update data counter <n> and add <dataLine> to <dataTime>
        n        = n + nLine;
        dataTime = [dataTime, dataLine];        
    end
    % break if <dataTime> is empty; error if <dataTime> has incorrect size
    if isempty(dataTime)
        break;
    elseif length(dataTime) < nPerTime
        error = 5;
        msg   = 'Error while reading data. Data has incorrect format.';
        return;
    end
    % update time counter <t>
    t  = t + 1;
    % update <time> if variable is time-dependent
    if timeVar
        time = [time,dataTime(1)];
    end
    % add single time step data to <data>
    data(t,1:nPerTime) = dataTime(timeVar+1:end);
    % reset counter and temporary data array <dataTime>
    n        = 0;
    dataTime = [];
end

%% POSTPROCESSING

% Reshape <data> to correct dimensions, in the right order
if timeVar
    data = reshape(data,[t,fliplr(dims)]);              % dimensions
    data = permute(data,[1,fliplr(2:numel(dims)+1)]);   % order
else
    data = reshape(data,fliplr(dims));                  % dimensions
    data = permute(data,fliplr(1:numel(dims)));         % order
end

% Close file
fclose(fid);

