function [c,ceq] = gen_con(z, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% alpha = param(1);
% beta = param(2);
% lambda_t = param(3);
% mx = param(4);
% N = param(5);
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;
mx = 6;
N = 50;
c = zeros(N,1);
for k = 1:N
    c(k) = alpha*exp(-beta*(z(1+(k-1)*mx) - lambda_t)^2) - z(5+(k-1)*mx);
end
ceq = [];
end

