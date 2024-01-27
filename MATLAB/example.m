clear all
close all
clc
format longEng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE USAGE OF MATLAB VERSION OF SDPEN PACKAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0     = [-0.5; 1.0];
lb     = [-0.5; -inf];
ub     = [ 0.5;  inf];
[f, c] = fob_con(x0);
n      = length(x0);
m      = length(c);
eps    = ones(m,1);

fprintf('------------------------\n');
fprintf('---- Initial values ----\n');
fprintf('------------------------\n');
fprintf(' fob = %e\n\n',f);

for i = 1:n
    fprintf(' x(%d) = %e\n',i,x0(i));
end
fprintf('\n');

for i = 1:m
    if (max(0.0,c(i)) < 1.0)
        eps(i) = 1.e-3;
    else
        eps(i) = 1.e-1;
    end
end

for i = 1:m
    fprintf(' con(%d) = %e   eps(%d) = %e\n',i,c(i),i,eps(i));
end
fprintf('------------------------\n');

epsiniz     = eps;
finiz       = f;
violiniz    = max(0.0,max(c));
alfa_stop   = 1.e-6;
nf_max      = 5000;
iprint      = -1;

fprintf('\n Now calling the optimizer ...\n');
[x, f, ni, nf, eps, istop] = sdpen(@fob_con,x0,lb,ub,alfa_stop,nf_max,iprint,eps);
fprintf(' ... done\n\n');

[fob, constr] = fob_con(x);

fprintf('\n------------------------\n');
fprintf('----  Final values -----\n');
fprintf('------------------------\n');
fprintf('\n');
fprintf('fob = %e\n',fob);
fprintf('\n');
for i = 1:n
    fprintf(' x(%d) = %e\n',i,x(i));
end
fprintf('\n');

for i = 1:m
    fprintf(' con(%d) = %e   eps(%d) = %e\n',i,constr(i),i,eps(i));
end
fprintf('------------------------\n');