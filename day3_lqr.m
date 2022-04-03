run('day2_optimalpath.m')

Q_lqr = diag([1 1 1 1]); %Four states
R_lqr = 0.01; %One input

K = dlqr(A_d,B_d,Q_lqr,R_lqr);

x_opt = timeseries([x1 x2 x3 x4]',t);