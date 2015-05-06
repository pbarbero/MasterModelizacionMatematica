function x = dtspg( T, b )
%function x = dtspg( T, b )

msg = nargchk(2, 2, nargin);
if size(msg) ~= 0
  error(msg); 
end

if ~isreal(T)
  error('T parameter in tpsymm must be a real matrix');
end

if ~isreal(b)
  error('b parameter in tpsymm must be a real vector');
end
if size(T,2) ~= size(b,1)
  error('the length of t and b parameters must be equal');
end
x = mex_dtspg(T, b);
