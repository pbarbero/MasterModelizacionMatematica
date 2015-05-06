% Common test data
n=10;
x=randn(n,1);
bSize = 2;
nBlocks = n/2;

% Symmetric Topelitz test:  T_sy x = b_sy
t=randn(n,1);
T_sy = toeplitz( t );
b_sy = T_sy * x;
x_dts = dts( t, b_sy );
disp(sprintf('Relative error in dts:   ||x-x_dts|| / ||x|| = %7.1E', ...
              norm(x-x_dts)/norm(x)));

% Non-symmetric Toeplitz test: T x = b
t0 = randn(1,n);
t1 = randn(1,n);
t1(1) = t0(1);
b_dt = toeplitz(t0,t1) * x;
x_dt = dt(t0,t1,b_dt);
disp(sprintf('Relative error in dt:    ||x-x_dt|| / ||x|| = %7.1E', ...
                norm(x-x_dt)/norm(x)));


% Symmetric tridiagonal Topelitz test: T_trid x = b_trid
t0 = randn(1,1);
t1 = randn(1,1);
T_trid = toeplitz( [ t0; t1; zeros(n-2,1) ] );
b_trid = T_trid*x;
x_dpis = dpis(t0, t1, b_trid);
disp(sprintf('Relative error in dpis:  ||x-x_dpis|| / ||x|| = %7.1E', ...
              norm(x-x_dpis)/norm(x)));


% Symmetric Block Toeplitz test: T_dpis x = b_dpis
H = zeros( (2*nBlocks-1)*bSize,nBlocks*bSize );
for i=1:nBlocks
        tmp = round(rand(bSize)*10);
        for j=1:nBlocks
                H((i+j-2)*bSize+1:(i+j-2)*bSize+bSize, (j-1)*bSize+1:(j-1)*bSize+bSize) = tmp;
        end
end
T_dtspg =H'*H;
b_dtspg = T_dtspg * x;
x_dtspg = dtspg(T_dtspg(1:bSize,:), b_dtspg);
disp(sprintf('Relative error in dtspg: ||x-x_dtspg|| / ||x|| = %7.1E', ...
		norm(x-x_dtspg)/norm(x)));
