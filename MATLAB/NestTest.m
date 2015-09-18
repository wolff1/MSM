function [et, rel_err] = NestTest(M)
%NESTTEST assume cubic b-splines
%   test creating fine grid from coarse grid
%   test solving for coarse grid from fine grid

    p = 4;
    J = [3./4.; 1./2.; 1./8.;];
    %C = (rand(M,1)-0.5).*2000;
    C = ones(M,1);
    %C = zeros(M,1);
    %C(floor(M/2)+1) = 1.0;
    %C = [1:M]';
    f = zeros(2*M+1+p/2,1);

    % Compute fine grid coefficients from coarse grid coefficients
    for i = 0:M-1
        for j = p/2:-1:1
            f(2*i-j+p/2+1) = f(2*i-j+p/2+1) + C(i+1)*J(j+1);
            %fprintf('(%d,%d) -> C(%d) = %f, J(%d) -> %f, f(%d) -> %f\n', i, -j, i, C(i+1), -j, J(j+1), 2*i-j, f(2*i-j+p/2+1));
        end
        for j = 0:p/2
            f(2*i+j+p/2+1) = f(2*i+j+p/2+1) + C(i+1)*J(j+1);
            %fprintf('(%d,%d) -> C(%d) = %f, J(%d) -> %f, f(%d) -> %f\n', i, j, i, C(i+1), j, J(j+1), 2*i+j, f(2*i+j+p/2+1));
        end
        %fprintf('\n');
    end
    %f

    t1 = tic();
    % Build linear system and solve
    Jp = [2.; 11./8.; 1./8.;];
    A = toeplitz([Jp' zeros(1, M-length(Jp))]);
    fp = zeros(M,1);
    for i = 0:M-1
        c = 2*i+p/2+1;
        for j = -p/2:p/2
            fp(i+1) = fp(i+1) + f(c+j);
        end
    end
    %fp

    %   Solve for coarse grid coefficients
    %   NOTE: This appears to be incorrect!
    C1 = A\fp;
    et = toc(t1);
    %C1
    %fp1 = A*C;
    %fp1
    %fprintf('|fp-fp1| = %e\n', norm(fp-fp1));

    % Fancy stuff (solve related linear system)
    c = 2.0*ones(1,length(Jp)-1)*Jp(2:end) - Jp(1) + 1.0;
    Ap = A + c*eye(M);
    L = myCholesky(Ap,0,length(Jp));

%     %y = zeros(M,1);
%     %y(1) = 1.0;
%     y = ones(M,1);
%     b1 = Ap*y;
%     bpp = rand(M,1);
%     bp = fp + bpp;
%     b0 = bp - b1;
% 
%     tmp = L\b0;
%     x = L'\tmp;
%     x
% 
%     tmp = L\bp;
%     x2 = L'\tmp - y;
%     x2
% 
%     norm(C-x)/norm(C)
%     norm(x-x2)/norm(C)

%     tmp = L\fp;
%     y = L'\tmp;
%     x = (eye(M)-c*inv(Ap))\y;
%     norm(C-x)/norm(C)

    Ap2 = Ap*Ap - c*c*eye(M);
    b = Ap\fp;
    Apb = Ap*Ap*b + c*Ap*b;
    x = Ap2\Apb;
    norm(C-x)/norm(C)

%     % Use recurrence to calculate coarse grid coefficients
%     C2 = zeros(M,1);
%     for i = 0:p/2-1 % index of current C value
%         % Calculate f'_{-p/2}
%         sump = 0.0;
%         sumn = 0.0;
%         for j = 1:2*i+1
%             if f(j) > 0.0
%                 sump = sump + f(j);
%             else
%                 sumn = sumn + f(j);
%             end
%         end
%         C2(i+1) = sump + sumn;
% 
%         sump = 0.0;
%         sumn = 0.0;
%         for j = i-1:-1:0    % j is # of terms in sum
%             if C2(j+1) >= 0.0
%                 sump = sump + Jp(abs(p/2-(i-j))+1)*C2(j+1);
%             else
%                 sumn = sumn + Jp(abs(p/2-(i-j))+1)*C2(j+1);
%             end
%         end
% 
%         if C2(i+1) >= 0.0
%             C2(i+1) = 8.0*((C2(i+1)-sumn)-sump);
%         else
%             C2(i+1) = 8.0*((C2(i+1)-sump)-sumn);
%         end
%     end
% 
% %     for i = p/2:M-1 % index of current C value
%     for i = p/2:M-2 % index of current C value
%         C2(i+1) = fp(i-p/2+1);
%         sump = 0.0;
%         sumn = 0.0;
%         for j = i-1:-1:max(i-p,0)
%             if C2(j+1) >= 0.0
%                 sump = sump + Jp(abs(p/2-(i-j))+1)*C2(j+1);
%             else
%                 sumn = sumn + Jp(abs(p/2-(i-j))+1)*C2(j+1);
%             end
%         end
%         if C2(i+1) >= 0.0
%             C2(i+1) = 8.0*((C2(i+1)-sumn)-sump);
%         else
%             C2(i+1) = 8.0*((C2(i+1)-sump)+sumn);
%         end
%     end
% 
%     % C2_{M-1} is special (i.e. need to "collapse" remaining equations in linear system into one)
%     i = M-1;
%     C2(i+1) = fp(i-p/2+1:end)'*ones(p/2+1,1); % added up to fp(i)
%     sump = 0.0;
%     sumn = 0.0;
%     for j = i-1:-1:max(i-p,0)
%         if C2(j+1) >= 0.0
%             coef = ones(1,p/2+1)*A(M-p/2:M,j+1);
%             sump = sump + (coef)*C2(j+1);
%         else
%             coef = ones(1,p/2+1)*A(M-p/2:M,j+1);
%             sumn = sumn + (coef)*C2(j+1);
%         end
%     end
%     if C2(i+1) >= 0.0
%         C2(i+1) = 2.0*((C2(i+1)-sumn)-sump)/7.0;
%     else
%         C2(i+1) = 2.0*((C2(i+1)-sump)+sumn)/7.0;
%     end
%     fprintf('|C-C1|/|C| = %e, |C-C2|/|C| = %e, |C1-C2|/|C| = %e\n', norm(C-C1)/norm(C), norm(C-C2)/norm(C), norm(C1-C2)/norm(C));
    rel_err = norm(C-C1)/norm(C);
    fprintf('|C-C1|/|C| = %e, et=%f\n', rel_err, et);
end
