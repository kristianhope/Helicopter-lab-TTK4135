run ('day2_optimalpath.m')

cd ('data_files/day2/round2')

load('day_2_r_01.mat');
z_real1 = state(:,:);
load('day_2_r_1.mat');
z_real2 = state(:,:);
load('day_2_r_10.mat');
z_real3 = state(:,:);

load('u_01.mat');
u_01 = u;
load('u_1.mat')
u_1 = u;
load('u_10.mat')
u_10 = u;

t_sim1 = z_real1(1,:);

extrainp_titles_labels = {'interpreter','latex','fontsize',22};
extrainp_legend = {'interpreter','latex','fontsize',15}; 

%% Plot originally from matlabscript day2_optimalpath.m
figure(2)
subplot(511)
stairs(t,u),grid
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(515)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')

%% Travel plot
figure();
hold on
grid on
%plot(t, x1); %%Optimal path
plot(t_sim1, z_real1(2,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(2,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(2,:)*pi/180,'LineWidth',1);
axis([-inf inf -inf inf])
set(gca, 'FontSize', 15);
legend({'$\lambda$, $q = 0.1$','$\lambda$, $q = 1$','$\lambda$, $q = 10$'},'Location','northeast',extrainp_legend{:});
%title("Travel w/ Open Loop Optimal Controller",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$\lambda$ [rad]",extrainp_titles_labels{:});


cd ('../../..')

%% Pitch plot

figure();
hold on
grid on
%plot(t, x1); %%Optimal path
plot(t_sim1, z_real1(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(4,:)*pi/180,'LineWidth',1);
axis([-inf inf -inf inf])
set(gca, 'FontSize', 15);
legend({'q = 0.1','q = 1','q = 10'},'Location','northeast',extrainp_legend{:});
%title("Pitch w/ Open Loop Optimal Controller",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$p$ [rad]",extrainp_titles_labels{:});

%% u plot

figure();
hold on
grid on
stairs(t, u_01,'LineWidth',1);
stairs(t, u_1,'LineWidth',1);
stairs(t, u_10,'LineWidth',1);
axis([-inf inf -inf inf])
set(gca, 'FontSize', 15);
legend({'$u$, $q = 0.1$','$u$, $q = 1$','$u$, $q = 10$'},'Location','northeast',extrainp_legend{:});
%title("Optimal input",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$u / p_c$ [rad]",extrainp_titles_labels{:});

%% comparison between opt and real

figure();
hold on
grid on
plot(t(1:100), x1(2:101),'LineWidth',1);
plot(t_sim1, z_real2(2,:)*pi/180,'LineWidth',1,'LineStyle','--');
axis([-inf inf -inf inf])
set(gca, 'FontSize', 15);
legend({'$\lambda*$','$\lambda$, $q = 1$'},'Location','northeast',extrainp_legend{:});
%title("Optimal vs. actual trajectory",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$\lambda$ [rad]",extrainp_titles_labels{:});

