
delta_t = 0.25;

%Discrete model
A_c = [0 1 0 0;
       0 0 -K_2 0;
       0 0 0 1;
       0 0 -K_1*K_pp -K_1*K_pp];
A1 = eye(4) + delta_t*A_c;
B_c = [0 0 0 K_1*K_pp]';
B1 = delta_t*B_c;

Q = diag([1 1 1 1]);

R = [1];

K = dlqr(A1,B1,Q,R);

x_opt = [x1; x2; x3; x4;];

u_opt = 1;