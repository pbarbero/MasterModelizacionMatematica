function x = dts( t, b, nb, piv )
%function x = dts( t, b[, nb[, piv]] )
%
% Solves the symmetric Toeplitz system T x = b
%
% Parameters:
%   t:   first column (or row) of the symmetric Toeplitz matrix T
%   b:   right hand side of the equation T x = b
%   nb:  block size of the algorithm (optional, default value is 64)
%   piv: use local diagonal pivoting if piv ~= 0 (optional, default value is 0) 
%   x:   solution of T x = b


msg = nargchk(2, 4, nargin);
if size(msg) ~= 0
  error(msg); 
end

if ~isreal(t)
  error('t parameter in dts must be a real vector');
end

if ~isreal(b)
  error('b parameter in dts must be a real vector');
end

if length(t) ~= length(b)
  error('the length of t and b parameters must be equal');
end

if nargin == 2
  nb = 64;
else 
  if ~isinteger(nb)
    error( 'nb parameter in dts must be an integer');
  end
  if ~(nb > 0) 
    error('block size in dts must be greater than 0');
  end 
end

if nargin <= 3
  piv = 0;
else
  if ~isintenger(piv)
    error('piv parameter in dts must be an integer');
  end
end

x = mex_dts(t, b, nb, piv);
