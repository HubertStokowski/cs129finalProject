function vals = getparam(fn,varargin)
vals = [];
for k = 1:length(varargin)
    vals(end+1) = sscanf(fn(findstr(fn,varargin{k}):end),[ varargin{k} '%f' ]);
end