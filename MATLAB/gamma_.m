% gamma_.m - This routine computes gamma(x) and gamma'(x)

% USAGE: gamma_(acc,obj,x) where:
%     acc is the accuracy parameters object
%     obj is the gamma parameters object
%     x is the independent variable

% DEPENDENCIES:   n/a

function [result,deriv] = gamma_(acc, obj, x)
% result = gamma(x), deriv = gamma'(x)
   if x >= 1.0
      result = 1.0/x;
      deriv = -1.0/x^2;
   else
      g0 = obj.g0;
      s = acc.s;
      k = acc.k_smooth;
      if (s > 1)
         % s>1: Break-point method
         bps = obj.bps;
         result = g0(1:k)'*x.^(2*(0:k-1)');
         deriv = (g0(2:k).*(2:2:2*k-2)')'*x.^((1:2:2*k-3)');
         deriv_tmp = 0.0;
         i = 1;
         while x >= bps(i)
            result = result + g0(i+k)*(x-bps(i))^(2*k);
            deriv_tmp = deriv_tmp + g0(i+k)*(x-bps(i))^(2*k-1);
            i = i + 1;
         end
         deriv = deriv + 2*k*deriv_tmp;
      else
         % s=1: Taylor method
         result = g0'*x.^((0:2:2*k)');
         deriv = (g0(2:end).*(2:2:2*k)')'*x.^((1:2:2*k-1)');
      end
   end

% End of file