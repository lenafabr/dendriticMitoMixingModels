function Sout = add2Struct(S,varargin)
% add the specified field names to the structure S

Sout = S;

for vc = 1:2:length(varargin)
    Sout.(varargin{vc}) = varargin{vc+1};
end   

end