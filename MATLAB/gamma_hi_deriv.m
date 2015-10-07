% gamma_hi_deriv.m - Given the objects which define gamma, this function
%                    returns the highest derivative for independent variable x

% USAGE: gamma_hi_deriv(acc, obj, x) where:
%     acc is the accuracy object
%     obj is the object defining gamma
%     x is a displacement

% DEPENDENCIES:   n/a

function [result] = gamma_hi_deriv(acc, obj, x)
   s = acc.s;
   k = acc.k_smooth;
   g0 = obj.g0;
   if x >= 1.0
      result = 1.0/x^(2*k+1);
   else
      result = 0.0;
      if s > 1
         % s>1: Break-point method
         i = 1;
         while x >= obj.bps(i)
            result = result + g0(i+k);
            i = i + 1;
         end
      else
         % s=1: Taylor Method
%         g0
%         obj.taylor_hi_deriv
         coef = g0 .* obj.taylor_hi_deriv;
%         obj.taylor_exp
         idx = obj.taylor_exp >= 0;
%         coef(idx)
%         obj.taylor_exp(idx)
         result = coef(idx)' * x.^obj.taylor_exp(idx);
%          result = result / factorial(k+1);
%          result
      end
   end
   
% End of file