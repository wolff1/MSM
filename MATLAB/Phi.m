% Phi.m  - Centered B-spline Phi(u) = Q_p(u + p) which has degree p-1
%        - This acts as parent basis function.

% USAGE: Phi(p,x) where
%           p determines the order of accuracy
%           x is the independent variable

% DEPENDENCIES:   N/A

function [fx, df_dx] = Phi(p, x)
   fx = 0;                 % Phi(x)
   df_dx = 0;              % Phi'(x)
   if abs(x) < p/2.0       % This can be removed if |x| is guaranteed to be < p/2
      Q = zeros(p,p);      % Q(u)
      Qp = zeros(p,p);     % Q'(u)
      u = x + p/2.0;       % b/c centered B-spline

      % start with k=1 (Q_1 is indicator function on [0,1])
      k = 1;
      for v = 0:p-k
          Q(1,v+1) = ((u-v) < 1.0 && (u-v) >= 0.0);
      end

      % use recurrence to build up to k=p-1
      for k = 2:p
          for v = 0:p-k
              Qp(k,v+1) = Q(k-1,v+1) - Q(k-1,v+2);
              Q(k,v+1) = (k*Q(k-1,v+2) + (u-v)*Qp(k,v+1))/(k-1);
          end
      end

      % Copy our computed value into the return values
      fx = Q(p,1);
      df_dx = Qp(p,1);
   end

function [fx, df_dx] = Phi2p(p, x)
   fx = 0;                 % Phi(x)
   df_dx = 0;              % Phi'(x)
   if abs(x) < p           % This can be removed if |x| is guaranteed to be < p
      Q = zeros(2*p,2*p);  % Q(u)
      Qp = zeros(2*p,2*p); % Q'(u)
      u = x + p;           % b/c centered B-spline

      % start with k=1 (Q_1 is indicator function on [0,1])
      k = 1;
      for v = 0:2*p-k
          Q(1,v+1) = ((u-v) < 1.0 && (u-v) >= 0.0);
      end

      % use recurrence to build up to k=2p-1
      for k = 2:2*p
          for v = 0:2*p-k
              Qp(k,v+1) = Q(k-1,v+1) - Q(k-1,v+2);
              Q(k,v+1) = (k*Q(k-1,v+2) + (u-v)*Qp(k,v+1))/(k-1);
          end
      end

      % Copy our computed value into the return values
      fx = Q(2*p,1);
      df_dx = Qp(2*p,1);
   end
   
% End of file
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         