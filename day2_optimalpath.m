%% Initialization and model definition
clear all
init08;
addpath('utilities');
% Discrete time system model. x = [lambda r p p_dot]'
dt	= 0.25; % sampling time
A_c = [0 1 0 0;
       0 0 -K_2 0;
       0 0 0 1;
       0 0 -K_1*K_pp -K_1*K_pd];
B_c = [0 0 0 K_1*K_pp]';

mx = size(A_c,2); % Number of states (number of columns in A)
mu = size(B_c,2); % Number of inputs(number of columns in B)

A_d = eye(mx) + dt*A_c;
B_d = dt*B_c;

x0 = [pi 0 0 0]';                       % Initial values

% Time horizon and initialization
N  = 100;                               % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Bounds
ul 	    = -pi/6;                        % Lower bound on control
uu 	    = pi/6;                         % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       =  gen_constraints(N,M,xl,xu,ul,uu);
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Generate the matrix G and the vector c (objecitve function weights in the QP problem) 
Q = diag([1 0 0 0]);
R = 0.1;                                           % Weight on input 1
I_N = eye(N);
G = 2*blkdiag(kron(I_N, Q), kron(I_N, R));        % Generate G

%% Generate system matrixes for linear model
Aeq = gen_aeq(A_d,B_d,N,mx,mu);         % Generate A
beq = [A_d*x0; zeros((N-1)*mx,1)];      % Generate b

%% Solve QP problem with linear model
tic
[z,lambda] = quadprog(G,[],[],[],Aeq,beq,vlb,vub,x0);
t1=toc;

%% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

padding_time = 5;                       % Seconds of padding
num_variables = padding_time/dt;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

t = 0:dt:dt*(length(u)-1);
u_t = timeseries(u,t);