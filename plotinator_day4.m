run ('day4_nonlinConstraint.m')

cd ('data_files/day4')
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

param = [alpha beta lambda_t mx length(z_real1(1,:))];
c_real1 = gen_con(z_real1, param);
c_real2 = gen_con(z_real2, param);
c_real3 = gen_con(z_real3, param);
c_real4 = gen_con(z_real4, param);
c_real5 = gen_con(z_real5, param);


t_sim1 = z_real1(1,:);

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
plot(t, c);
hold on
plot(t_sim1, z_real1(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real2(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real3(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real4(6,:)*pi/180,'LineWidth',1);
plot(t_sim1, z_real5(6,:)*pi/180,'LineWidth',1);
axis([0 20 -inf inf])
set(gca, 'FontSize', 15);
legend({'Constraint','Test 1','Test 2','Test 3','Test 4','Test 5'},'Location','northeast',extrainp_legend{:});
title("Elevation w/ Constraint",extrainp_titles_labels{:});
xlabel("$t$ [s]",extrainp_titles_labels{:});
ylabel("$e$ [rad]",extrainp_titles_labels{:});


cd ('../..')