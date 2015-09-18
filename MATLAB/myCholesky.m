% myCholesky.m - A suite of Cholesky factorizations.

function [L,Lp] = myCholesky(A,v,p)
%    t1 = tic();
%    [~,L1] = generalCholesky(A);
%    fprintf('t1 = %f\n', toc(t1));
% 
%    t2 = tic();
%    [~,L2] = bCholesky(A,p);
%    fprintf('t2 = %f\n', toc(t2));
% 
%    t3 = tic();
    L = btCholesky(A,v,p);
    Lp = 0;
%    fprintf('t3 = %f\n', toc(t3));
% 
%    t4 = tic();
%    L4 = ibtCholesky(v,p);
%    fprintf('t4 = %f\n', toc(t4));

%   t5 = tic();
%   [L,Lp] = sibtCholesky(v,p);
%   fprintf('t5 = %f\n', toc(t5));

%    % Display difference
%    fprintf('\tDIFF1: %e\n', norm(L3 - L1) / norm(L1));
%    fprintf('\tDIFF2: %e\n', norm(L3 - L2) / norm(L2));
   dummy = 0;

function [L,Lp] = sibtCholesky(v,p)
   % compute Cholesky factorization of matrix A
   % input:  v is row vector containing the values within the band, 1=diagonal
   %         p is bandwidth
   %         A is an singly infinite, symmetric, banded toeplitz matrix
   %           which has been folded twice and divided by 2 to exploit
   %           symmetry within the system
   % output: L is lower triangular matrix such that L * L' = A
   %           (represented as the vector containing nonzero elements in
   %           bandwidth)

   % Call ibtCholesky to get result of infinite Cholesky factorization
   L = ibtCholesky(v,p);

   % Calculate Cholesky for lower-right portion of coefficient matrix
   v = v(1:p);
   bw = [v(p:-1:1)' v(2:p)']; % full bandwidth vector
   A = zeros(2*p-1,p); % FAKE COEFFICIENT MATRIX
   A(1:p,:) = ones(p,1)*v(p:-1:1)';
   stencil = ones(p,1)*L;
   for i = 1:p-1
      % tmp is bw with i-most elements folded over
      tmp = bw(1:2*p-1-i) + [zeros(1,2*(p-i-1)) bw(1:i) 0.0];
      tmp = tmp(1:p);
      if (i == p-1), tmp = tmp / 2.0; end
      A(p+i,:) = tmp;

      % Begin Cholesky algorithm
      Lii_sum = 0.0;
      for j = 1:p-1
         Lij_sum = 0.0;
         for k = 1:j-1
            Lij_sum = Lij_sum + stencil(p,k)*stencil(j,k-j+p);
         end
         stencil(p,j) = (A(p+i,j) - Lij_sum) / stencil(j,p);% A'(p+i,j) = A(i,j)
         Lii_sum = Lii_sum + stencil(p,j)*stencil(p,j);
      end
      stencil(p,p) = sqrt(A(p+i,p) - Lii_sum);              % A'(p+i,p) = A(i,i)
      stencil(1:p-1,:) = stencil(2:p,:);                    % Roll over
   end
   
   % Return last p-1 rows of singly infinite coefficient matrix
   Lp = stencil(1:p-1,:);

function L = ibtCholesky(v,p)
   % compute Cholesky factorization of matrix A
   % input:  A is an infinite, symmetric, banded toeplitz matrix
   %         v is row vector containing the values within the band, 1=diagonal
   %         p is bandwidth
   % output: L is lower triangular matrix such that L * L' = A
   %           (represented as the vector containing nonzero elements in
   %           bandwidth)

   stencil = zeros(p,p);
   v = v(1:p);

   gmm1 = zeros(1,p);
   gm = zeros(1,p);
   gm_save_flag = 0;
   gm_save = zeros(1,p);
   last_diff = 0.0;

   % p steps of symmetric, banded toeplitz Cholesky to build helper table
   L = zeros(p,p);
   L(1,1) = sqrt(v(1));                               % v(1) = A(1,1)
   stencil(1,p) = sqrt(v(1));                          % v(1) = A(1,1)
   for i = 2:p-1
      imin = max(1,i-p+1);
      Lii_sum = 0.0;
      for j = imin:i-1
         Lij_sum = 0.0;
         for k = imin:j-1
            Lij_sum = Lij_sum + L(i,k)*L(j,k);
         end
         L(i,j) = (v(i-j+1) - Lij_sum) / L(j,j);      % v(i-j+1) = A(i,j)
         Lii_sum = Lii_sum + L(i,j)*L(i,j);
      end
      L(i,i) = sqrt(v(1) - Lii_sum);                  % v(1) = A(i,i)
      stencil(i,p+imin-i:p) = L(i,imin:i);
   end

   % Now run infinite Cholesky once helper table has been tabulated
   for i = 1:1000
      Lii_sum = 0.0;
      for j = 1:p-1
         Lij_sum = 0.0;
         for k = 1:j-1
            Lij_sum = Lij_sum + stencil(p,k)*stencil(j,k-j+p);
         end
         stencil(p,j) = (v(p-j+1) - Lij_sum) / stencil(j,p);% v(p-j+1) = A(i,j)
         Lii_sum = Lii_sum + stencil(p,j)*stencil(p,j);
      end
      stencil(p,p) = sqrt(v(1) - Lii_sum);                  % v(1) = A(i,i)

      % TOLERANCE CHECK
      if gm_save_flag
         this_diff = norm(gm_save - stencil(p,:));
         if this_diff <= last_diff
            %STOP
%            fprintf('BREAKING B/C THIS_DIFF <= LAST_DIFF, ITERATION<%d>\n',i);
            break;
         else
            % SET LAST DIFF
%            fprintf('SETTING LAST_DIFF = THIS_DIFF\n');
            last_diff = this_diff;
         end
      else
         gmm1 = gm;
         gm = stencil(p-1,:);
         gmp1 = stencil(p,:);
         norm_val1 = norm(gmp1 - gm);
         norm_val2 = norm(gm - gmm1);
         if (norm(gmp1 - gm) >= norm(gm - gmm1))
            % SAVE OFF GM
%            fprintf('SAVING OFF G_M\n');
            gm_save = gm;
            gm_save_flag = 1;
         else
            % KEEP GOING
%            fprintf('FIRST CONVERGENCE CRITERIA NOT MET. <%e> <%e>\n', ...
%                                                     norm_val1, norm_val2);
         end
      end

      % Roll over
      stencil(1:p-1,:) = stencil(2:p,:);
   end

   % Return last row of stencil as it has the converged values for bandwidth
   L = stencil(p,:);

function L = btCholesky(A,v,p)
   % compute Cholesky factorization of matrix A
   % input:  A is symmetric, banded toeplitz matrix
   %         v is row vector containing the values within the band, 1=diagonal
   %         p is bandwidth
   % output: L is lower triangular matrix such ath L * L' = A

   n = size(A,1);
   L = zeros(n,n);

   L(1,1) = sqrt(A(1,1));
   for i = 2:n
      imin = max(1,i-p+1);
      Lii_sum = 0.0;
      for j = imin:i-1
         Lij_sum = 0.0;
         for k = imin:j-1
            Lij_sum = Lij_sum + L(i,k)*L(j,k);
         end
         L(i,j) = (A(i,j) - Lij_sum) / L(j,j);
         Lii_sum = Lii_sum + L(i,j)*L(i,j);
      end
      L(i,i) = sqrt(A(i,i) - Lii_sum);
   end

function [s, L] = bCholesky(A,p)
   % compute Cholesky factorization of matrix A
   % input:  A is banded & symmetric matrix
   %         p is bandwidth
   % output: L is lower triangular matrix such ath L * L' = A
   %         s =  1 if A is s.p.d
   %             -1 if A is s.n.d
   %              0 if A is neither s.p.d nor s.n.d

   n = size(A,1);
   s = sign(A(1,1));
   L = zeros(n,n);

   for i = 1:n

       % compute the diagonal element of column i
       sum = 0;
       for j = max(1,i-p+1):i-1
           sum = sum + L(i,j)*L(i,j);
       end
       if(s * A(i,i) - sum <= 0) % neither case
           s = 0;
           fprintf('input matrix is neither s.p.d nor n.p.d\n');
           return
       end
       L(i,i) = sqrt(s * A(i,i) - sum);

       % compute elements below the diagonal of column i
       for j = i+1:min(i+p-1,n)
           sum = 0;
           for k = max(1,i-p+1):i-1
               sum = sum + L(i,k)*L(j,k);
           end
           L(j,i) = s*(A(i,j) - sum)/L(i,i);
       end
   end

function [s, L] = generalCholesky(A)
% MAW - I STOLE THIS FROM A CS515 SAMPLE SOLUTION
   % compute Cholesky factorization of matrix A
   % input:  A is symmetric matrix
   % output: L is lower triangular matrix such ath L * L' = A
   %         s =  1 if A is s.p.d
   %             -1 if A is s.n.d
   %              0 if A is neither s.p.d nor s.n.d

   n = size(A,1);
   s = sign(A(1,1));
   L = zeros(n,n);

   for i = 1:1:n

       % compute the diagonal element of column i
       sum = 0;
       for j = 1:1:i-1
           sum = sum + L(i,j)*L(i,j);
       end
       if(s * A(i,i) - sum <= 0) % neither case
           s = 0;
           fprintf('input matrix is neither s.p.d nor n.p.d\n');
           return
       end
       L(i,i) = sqrt(s * A(i,i) - sum);

       % compute elements below the diagonal of column i
       for j = i+1:1:n
           sum = 0;
           for k = 1:1:i-1
               sum = sum + L(i,k)*L(j,k);
           end
           L(j,i) = s*(A(i,j) - sum)/L(i,i);
       end
   end

% End of file
