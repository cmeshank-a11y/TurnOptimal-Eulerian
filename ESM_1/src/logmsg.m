function logmsg(varargin)
%LOGMSG  Print to the console and append to the run log.
%   logmsg('#open', file) opens the log, logmsg('#close') closes it,
%   logmsg(fmt, args...) prints like fprintf.
persistent fid
if nargin >= 1 && ischar(varargin{1}) && strcmp(varargin{1}, '#open')
    fid = fopen(varargin{2}, 'w');
    return
end
if nargin >= 1 && ischar(varargin{1}) && strcmp(varargin{1}, '#close')
    if ~isempty(fid) && fid > 0, fclose(fid); end
    fid = [];
    return
end
fprintf(varargin{:});
if ~isempty(fid) && fid > 0
    fprintf(fid, varargin{:});
end
end
