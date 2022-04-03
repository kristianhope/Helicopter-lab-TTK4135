%% Initialization and model definition
clear all
init08; % Change this to the init file corresponding to your helicopter
addpath('utilities');
% Discrete time system model. x = [lambda r p p_dot e e_dot]'
dt	= 0.25; % sampling time
A_c = [0 1 0 0 0 0;
       0 0 -K_2 0 0 0;
       0 0 0 1 0 0;
       0 0 -K_1*K_pp -K_1*K_pd 0 0;
       0 0 0 0 0 1;
       0 0 0 0 -K_3*K_ep -K_3*K_ed];
B_c = [0 0;
       0 0;
       0 0;
       K_1*K_pp 0;
       0 0;
       0 K_3*K_ep];

mx = size(A_c,2); % Number of states (number of columns in A)
mu = size(B_c,2); % Number of inputs(number of columns in B)

A_d = eye(mx) + dt*A_c;
B_d = dt*B_c;

x0 = [pi 0 0 0 0 0]';                   % Initial values

% Time horizon and initialization
N  = 40;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds
ul 	    = [-pi/6; -inf];                % Lower bound on control
uu 	    = [pi/6; inf];                  % Upper bound on control
xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                        % Lower bound on state x3 (pitch)
xu(3)   = uu(1);                        % Upper bound on state x3 (pitch)

% Generate constraints on measurements and inputs
[vlb,vub]       =  gen_constraints(N,M,xl,xu,ul,uu); 
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix G and the vector c (objecitve function weights in the QP problem) 
Q1 = 2*diag([1 0 0 0 0 0]);
P1 = zeros(mu,mu);
P1(1,1) = 10;                           % Weight on input 1
P1(2,2) = 10;                           % Weight on input 2
I_N = eye(N);
G = blkdiag(kron(I_N, Q1), kron(I_N, P1));        % Generate Q

%% Generate system matrixes for linear model
Aeq = gen_aeq(A_d,B_d,N,mx,mu);             % Generate A
beq = [A_d*x0; zeros((N-1)*mx,1)];          % Generate b

%% Inequality constraint for elevation
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;
param = [alpha beta lambda_t mx N];
nonlincon = @(z)gen_con(z,param);
%% Objective function

objfun = @(z) 0.5*z'*G*z;

% Use fmincon because of nonlinear constraint
opt = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',40000);
tic
[z, fval,exitflag,output] = fmincon(objfun,z0,[],[],Aeq,beq,vlb,vub,nonlincon,opt);
toc
c = gen_con(z,param) + z(5:mx:N*mx);
c = [c; 0;];

%% Extract control inputs and states
%%u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
u  = [z(N*mx+1:mu:end)'; z(N*mx+2:mu:end)'];
u = [u z(N*mx + M*mu -1:N*mx + M*mu)];
u1 = u(1,:);
u2 = u(2,:);

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

padding_time = 8;                       % Seconds of zero padding
num_variables = padding_time/dt;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1'; zero_padding];
u2   = [zero_padding; u2'; zero_padding];
u    = [u1'; u2'];

x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];
c = [zero_padding; c; zero_padding];

t = 0:dt:dt*(length(u1)-1);
u_t = timeseries(u,t);

%% LQR
Q_lqr = diag([10 15 10 10 10 15]); %four states
R_lqr = diag([0.1 0.05]); %two input

K = dlqr(A_d,B_d,Q_lqr,R_lqr);

x_optimal = [x1 x2 x3 x4 x5 x6]';

x_opt = timeseries(x_optimal,t);

