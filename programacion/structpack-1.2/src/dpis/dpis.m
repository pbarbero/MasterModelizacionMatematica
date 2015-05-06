function  x = dpis(t0, t1, b, method) 
% function  x = dpis(t0, t1, b[, method]) 
%
% Solves the real symmetric tridiagonal Topelitz system T x = b, where
%
%       /                       \
%      |  t0 t1 0  ...  0  0  0  |
%      |  t1 t0 t1 ...  0  0  0  |
%      |  0  t1 t0 ...  0  0  0  |
%  T = .  .  .  .  ...  .  .  .  .
%      |  0  0  0  ...  t0 t1 0  |
%      |  0  0  0  ...  t1 t0 t1 |
%      |  0  0  0  ...  0  t1 t0 |
%       \                       /
%
% The problem dimension is the dimension of vector b.
%
% Parameters:
%   t0: diagonal element of the symmetric tridiagonal matrix T
%   t1: super-(sub-)diagonal element of the symmetric tridiagonal matrix T
%   b: right hand side of the equation   T x = b
%   method (optional):  
%     - 'rojo': Use the Modified Rojo Method (default)
%     - 'dst':  Use the DST Method
%     - 'ldlt': Use the LDLt Method 
%


msg = nargchk(3, 4, nargin);
if size(msg) ~= 0
  error(msg);
end

if ~isreal(t0) || ~isscalar(t0)
  error('t0 parameter in dpis must be a real scalar');
end

if ~isreal(t1) || ~isscalar(t1)
  error(' t1 parameter in dpis must be a real scalar ');
end

if ~isreal(b) || ~isvector(b)
  error(' b parameter in dpis must be a real vector ');
end

if nargin == 3
  method = 'rojo'; 
end

[r, c] = size(b);
n=r*c;

if strcmp(method, 'rojo') || strcmp(method, 'dst') || strcmp(method,' ldlt')
  x = mex_dpis(n, t0, t1, b, method);
else
  error(sprintf("method %s not recognised", method));
end
