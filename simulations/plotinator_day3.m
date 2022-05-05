clear all
run ('day3_lqr.m')

cd ('data_files/day3')
load('day_3_test1.mat');
z_real1 = state(:,:);
load('day_3_test2.mat');
z_real2 = state(:,:);
load('day_3_test3.mat');
z_real3 = state(:,:);
load('day_3_test4.mat');
z_real4 = state(:,:);
load('day_3_test5.mat');
z_real5 = state(:,:);
load('day_3_test6.mat');
z_real6 = state(:,:);
load('day_3_test7.mat');
z_real7 = state(:,:);
load('day_3_test8.mat');
z_real8 = state(:,:);
load('day_3_test9.mat');
z_real9 = state(:,:);
load('day_3_test10.mat');
z_real10 = state(:,:);

t_sim1 = z_real1(1,:);

extrainp_titles_labels = {'interpreter','latex','fontsize',22};
extrainp_legend = {'interpreter','latex','fontsize',15}; 

%% Travel plot
figure();
hold on
plot(t, x1); %%Optimal path
plot(t_sim1, z_real10(2,:)*pi/180,'LineWidth',1);
axis([-inf 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'$\lambda^{*}$','$\lambda$'},'Location','northeast',extrainp_legend{:});
title("Travel w/ Open Loop Optimal Controller and LQR",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$\lambda$ [rad]",extrainp_titles_labels{:});

%% Pitch plot
figure();
hold on
plot(t, x3); %%Optimal path
plot(t_sim1, z_real10(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real1(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real4(4,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real5(4,:)*pi/180,'LineWidth',1);
axis([-inf 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'$p_c$','$p$'},'Location','northeast',extrainp_legend{:});
title("Pitch w/ Open Loop Optimal Controller and LQR",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$p$ [rad]",extrainp_titles_labels{:});

%%
cd ('../..')