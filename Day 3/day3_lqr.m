run('day2_optimalpath.m')

Q_lqr = diag([0.5 1 1 1]);
R_lqr = 1;
K = dlqr(A_d,B_d,Q_lqr,R_lqr);

x_opt = timeseries([x1 x2 x3 x4]',t);