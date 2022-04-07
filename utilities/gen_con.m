function [c,ceq] = gen_con(z, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
alpha = param(1);
beta = param(2);
lambda_t = param(3);
mx = param(4);
N = param(5);
c = zeros(N,1);
for k = 1:N
    c(k) = alpha*exp(-beta*(z((k-1)*mx+1) - lambda_t)^2) - z((k-1)*mx + 5);
end
ceq = [];
end

