function [ output_args ] = FindGridSpacing(tau, D, N)
%FINDGRIDSPACING Summary of this function goes here
%   Detailed explanation goes here
    output_args = 0;

    a = 1.0;        % FIXME - UNACCOUNTED FOR
    p = 4;
    while 1
        if p > 10
            break;
        end
        acc.s = 1;
        acc.k_smooth = p/2 - 1;    % for building gamma
        gamma_obj = gamma_init(acc);
        hi_deriv = gamma_hi_deriv(acc, gamma_obj, 0.0)/a;   % assumption that max high derivative is approximately the same as at x=0.0

        lead_term = 1.0;
        for i = 1:2:p-1
            lead_term = (lead_term * i)/(i+1);
        end

        fprintf('p = %d, 2gamma_n(p) - gamma_n(p/2) = %f, gamma_n'' = %f\n', p, 2.0*gamma_n(p)-gamma_n(p/2), gamma_n_p(p));
        p = p + 2;
        break;
    end

end

function sum = gamma_n(n)
    nat = 1:n;
    inv = 1.0./nat;
    sum = ones(1,n)*inv' - log(n);
end

function sum = gamma_n_p(n)
    one = 1.0./(1:(n/2));
    two = 1.0./((n/2+1):n);
    sum = ones(1,n/2)*one' + 2.0*(ones(1,n/2)*two') - log(2.0*n);
end