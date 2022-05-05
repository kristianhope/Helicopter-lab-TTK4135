clear all
run ('day4_nonlinConstraint.m')
cd ('data_files/day4/round2/')
load('day_4_test1.mat');
z_real1 = state(:,:);
load('day_4_test2.mat');
z_real2 = state(:,:);
load('day_4_test3.mat');
z_real3 = state(:,:);
load('day_4_test4.mat');
z_real4 = state(:,:);
load('day_4_test5.mat');
z_real5 = state(:,:);
load('day_4_test6.mat');
z_real6 = state(:,:);

cd ('../../..')
cd ('data_files/day4/round1/')
load("data_files/day4/round1/day_4_test1.mat")
z_real21 = state(:,:);
load('day_4_test2.mat');
z_real22 = state(:,:);
load('day_4_test3.mat');
z_real23 = state(:,:);
load('day_4_test4.mat');
z_real24 = state(:,:);
load('day_4_test5.mat');
z_real25 = state(:,:);

param = [alpha beta lambda_t mx length(z_real4(1,:))];
c_real1 = gen_con(z_real1, param);
c_real2 = gen_con(z_real2, param);
c_real3 = gen_con(z_real3, param);
c_real4 = gen_con(z_real4, param);
%c_real5 = gen_con(z_real5, param);


t_sim1 = z_real4(1,:);
t_sim2 = z_real21(1,:);

extrainp_titles_labels = {'interpreter','latex','fontsize',22};
extrainp_legend = {'interpreter','latex','fontsize',15}; 

%% Plot originally from matlabscript day4_lqr.m
figure(2)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5,'mo'),grid
hold on
plot(t,c,'m',t,c,'b'),grid
ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')

%% Elevation plot
figure();
% subplot(211);
% plot(t(1:115-20+1), c(20:end));
hold on
grid on
t_sim3 = 0:0.002:19.002;
% plot(t_sim1(1:12501-2500+1), z_real1(6,2500:end)*pi/180,'LineWidth',1);
% plot(t_sim1(1:12501-2500+1), z_real2(6,2500:end)*pi/180,'LineWidth',1);
% plot(t_sim1(1:12501-2500+1), z_real3(6,2500:end)*pi/180,'LineWidth',1);
plot(t_sim3, z_real1(6,2500:12001)*pi/180,'LineWidth',1);
%plot(t_sim1, z_real2(6,:)*pi/180,'LineWidth',1);
%plot(t_sim1, z_real3(6,:)*pi/180,'LineWidth',1);
plot(t_sim3, z_real6(6,2500:12001)*pi/180,'LineWidth',1);
plot(t_sim3, z_real3(6,2500:12001)*pi/180,'LineWidth',1);
plot(t_sim3, z_real2(6,2500:12001)*pi/180,'LineWidth',1);
plot(t(1:96-20+1), c(20:96),'LineWidth',1.3);

axis([0 19 -inf inf])
set(gca, 'FontSize', 15);
legend({'$e$, $\mathbf{Q}$ = diag([1 1 1 1 1 1])', ...
    '$e$, $\mathbf{Q}$ = diag([15 1 1 1 30 5])', ...
    '$e$, $\mathbf{Q}$ = diag([1 1 1 1 $10^3$ 1])', ...
    '$e$, $\mathbf{Q}$ = diag([1 1 1 1 $10^4$ 1])', ...
    'Constraint'},'Location','northeast',extrainp_legend{:});
%title("Elevation w/ Constraint",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$e$ [rad]",extrainp_titles_labels{:});
% subplot(212);
% hold on
% grid on
% % plot(t_sim3, z_real1(2,2500:10001)*pi/180,'LineWidth',1);
% % plot(t_sim3, z_real4(2,2500:10001)*pi/180,'LineWidth',1);
% % plot(t_sim3, z_real21(2,2500:10001)*pi/180,'LineWidth',1);
% plot(t_sim1, z_real1(2,1:end)*pi/180,'LineWidth',1);
% plot(t_sim1, z_real2(2,:)*pi/180,'LineWidth',1);
% plot(t_sim2, z_real21(2,:)*pi/180,'LineWidth',1);
% plot(t,x1);
% axis([0 20 -inf inf])
% set(gca, 'FontSize', 15);
% legend({'Test 1','Test 2','Test 3'},'Location','northeast',extrainp_legend{:});
% %title("Elevation w/ Constraint",extrainp_titles_labels{:});
% xlabel("$t$ [s]",extrainp_titles_labels{:});
% ylabel("$e$ [rad]",extrainp_titles_labels{:});
%% elev & travel & pitch
figure();
grid on
t_sim32 = 0:0.002:15.5;
subplot(311);
plot(t_sim32, z_real6(6,2500:10250)*pi/180,'LineWidth',1);
hold on
grid on
plot(t(1:83-21+1), x5(21:83),'LineWidth',1.3,'LineStyle','--');
axis([0 15.5 -inf inf])
set(gca, 'FontSize', 15);
legend({'$e$, $\mathbf{Q}$ = diag([15 1 1 1 30 5])','$e^*$'},'Location','northeast',extrainp_legend{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$e$ [rad]",extrainp_titles_labels{:});

subplot(312);
%plot(t_sim32, z_real6(2,2500:10250)*pi/180,'LineWidth',1);
plot(t_sim32,z_real6(2,2500:10250)*pi/180,'LineWidth',1);

hold on;
grid on
plot(t(1:83-21+1), x1(21:83),'LineWidth',1.3,'LineStyle','--');
%plot(t,x1);
axis([0 15.5 -inf inf])
set(gca, 'FontSize', 15);
legend({'$\lambda$, $\mathbf{Q}$ = diag([15 1 1 1 30 5])','$\lambda^*$'},'Location','northeast',extrainp_legend{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$\lambda$ [rad]",extrainp_titles_labels{:});

subplot(313);
%plot(t_sim32, z_real6(2,2500:10250)*pi/180,'LineWidth',1);
plot(t_sim32,z_real6(4,2500:10250)*pi/180,'LineWidth',1);
hold on;
grid on
plot(t(1:83-21+1), x3(21:83),'LineWidth',1.3,'LineStyle','--');
%plot(t,x3);
axis([0 15.5 -inf inf])
set(gca, 'FontSize', 15);
legend({'$p$, $\mathbf{Q}$ = diag([15 1 1 1 30 5])','$p^*$'},'Location','northeast',extrainp_legend{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$p$ [rad]",extrainp_titles_labels{:});

%% 2
figure();
%plot(t(1:115-20+1), c(20:end));
hold on
grid on
plot(t_sim1, z_real1(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(6,:)*pi/180,'LineWidth',1);
%plot(t_sim2, z_real23(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(6,:)*pi/180,'LineWidth',1);
axis([0 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'Constraint','Test 1','Test 2','Test 3','Test 4','Test 5'},'Location','northeast',extrainp_legend{:});
title("Elevation w/ Constraint",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$e$ [rad]",extrainp_titles_labels{:});


%% plot u_k

z_real1_chopped = z_real1(2:7, 1:125:end)*pi/180;
z_real2_chopped = z_real2(2:7, 1:125:end)*pi/180;
z_real3_chopped = z_real3(2:7, 1:125:end)*pi/180;

x_all = ([x1 x2 x3 x4 x5 x6]');

u_LQR1 = [];
u_LQR2 = [];
u_LQR3 = [];
for k = 1:length(z_real1_chopped)
    u_LQR1 = [u_LQR1, (u(:,k)-K*(z_real1_chopped(:,k)-x_all(:,k)))];
    u_LQR2 = [u_LQR2, (u(:,k)-K*(z_real2_chopped(:,k)-x_all(:,k)))];
    u_LQR3 = [u_LQR3, (u(:,k)-K*(z_real3_chopped(:,k)-x_all(:,k)))];
end
figure();
hold on
grid on
plot(t(1:101),u1(1:101))
plot(t(1:101),u_LQR1(1,:));
plot(t(1:101),u_LQR2(1,:));
plot(t(1:101),u_LQR3(1,:));

axis([0 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'ec','Test 1','Test 2','Test 3'},'Location','northeast',extrainp_legend{:});
title("Elevation w/ Constraint",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$p_c$ [rad]",extrainp_titles_labels{:});

figure();
hold on
grid on
plot(t(1:101),u2(1:101))
plot(t(1:101),u_LQR1(2,:));
plot(t(1:101),u_LQR2(2,:));
plot(t(1:101),u_LQR3(2,:));

axis([0 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'ec','Test 1','Test 2','Test 3'},'Location','northeast',extrainp_legend{:});
title("Elevation w/ Constraint",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$e_c$ [rad]",extrainp_titles_labels{:});


%% Pitch plot
figure();
%plot(t, x3);
hold on
plot(t_sim1, z_real1(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(4,:)*pi/180,'LineWidth',1);
%plot(t_sim1, z_real4(4,:)*pi/180,'LineWidth',1);
%plot(t_sim1, z_real5(4,:)*pi/180,'LineWidth',1);
%plot(t_sim1, z_real6(4,:)*pi/180,'LineWidth',1);
axis([0 20 -inf inf])
set(gca, 'FontSize', 15);
legend({},'Location','northeast',extrainp_legend{:});
title("Pitch",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$p$ [rad]",extrainp_titles_labels{:});
cd ('../../..')