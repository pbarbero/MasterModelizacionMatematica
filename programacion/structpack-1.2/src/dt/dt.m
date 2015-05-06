function x = dt(u, v, b, nb, piv )
%function x = dt(u, v, b[,nb [, piv]] )
%
% Solves the Toeplitz system T x = b
%
% Parameters:



msg = nargchk( 3, 5, nargin );
if size( msg )~=0
  error( msg ); 
end

if ~isreal( u )
  error( 'u parameter in dt must be a real vector' );
end

if ~isreal( v )
  error( 'v parameter in dt must be a real vector' );
end

if ~isreal( b )
  error( 'b parameter in dt must be a real vector' );
end

if length( u ) ~= length( b )
  error( 'the length of u and b parameters must be equal' );
end
if length( v ) ~= length( b )
  error( 'the length of v and b parameters must be equal' );
end

if ( nargin == 3 )
  nb = 64; 
else 
  if ~isinteger( nb )
    error( 'nb parameter in dt must be an integer' );
  end
  if ~( nb > 0 ) 
    error( 'block size in dt must be greater than 0' );
  end 
end

if ( nargin <= 4 )
  piv = 0;
else
  if ~isintenger( piv )
    error( 'piv parameter in dt must be an integer' );
  end
end

x = mex_dt( u, v, b, nb, piv );
