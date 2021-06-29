%%
clear
close all
clc

%% LTI systems and given parameters

LTI1.A=[1 0 2 0; 0 1 0 2; 0 0 3 0; 0 0 0 3];
LTI1.B=[2 0;0 2;3 0;0 3];
LTI1.x0=[-10;10;-1;1];

LTI2.A=[1 0 3 0; 0 1 0 3; 0 0 7 0; 0 0 0 7];
LTI2.B=[3 0; 0 3; 7 0; 0 7];
LTI2.x0=[10;10;1;1];

LTI3.A=[1 0 1 0; 0 1 0 1; 0 0 1.1 0; 0 0 0 1.1];
LTI3.B=[1 0; 0 1; 1.1 0; 0 1.1];
LTI3.x0=[10;-10;1;-1];

LTI4.A=[1 0 6 0; 0 1 0 6; 0 0 20 0; 0 0 0 20];
LTI4.B=[6 0;0 6;20 0; 0 20];
LTI4.x0=[-10;-10;-1;-1];

Tfinal=5;
umax=100;
N = 4; 

% Definition of system dimension
dim.nx=4;     %state dimension
dim.nu=2;     %input dimension
dim.N=Tfinal; %horizon

%Definition of quadratic cost function
weight.Q=eye(dim.nx);   %weight on output
weight.R=eye(dim.nu);   %weight on input

% Generation of prediction model 1
predmod1=predmodgen(LTI1,dim);            
[H1,h1]=costgen(predmod1,dim,LTI1);

% Generation of prediction model 2
predmod2=predmodgen(LTI2,dim);            
[H2,h2]=costgen(predmod2,dim,LTI2);

% Generation of prediction model 3
predmod3=predmodgen(LTI3,dim);            
[H3,h3]=costgen(predmod3,dim,LTI3);

% Generation of prediction model 4
predmod4=predmodgen(LTI4,dim);            
[H4,h4]=costgen(predmod4,dim,LTI4);

%% Constraints

%Equality constraints
%Constraints of model 1 
b_eq_1 = predmod1.T(dim.nx*(dim.N-1)+1:end,:)*LTI1.x0;
A_eq_1 = predmod1.S(dim.nx*(dim.N-1)+1:end,:);
%Constraints of model 2
b_eq_2 = predmod2.T(dim.nx*(dim.N-1)+1:end,:)*LTI2.x0;
A_eq_2 = predmod2.S(dim.nx*(dim.N-1)+1:end,:);
%Constraints of model 3
b_eq_3 = predmod3.T(dim.nx*(dim.N-1)+1:end,:)*LTI3.x0;
A_eq_3 = predmod3.S(dim.nx*(dim.N-1)+1:end,:);
%Constraints of model 4
b_eq_4 = predmod4.T(dim.nx*(dim.N-1)+1:end,:)*LTI4.x0;
A_eq_4 = predmod4.S(dim.nx*(dim.N-1)+1:end,:);
%Equality constraints (They are the same for every agent)
b_ineq = ones(2*dim.N*dim.nu,1)*umax/Tfinal;
A_ineq = blkdiag([1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1],[1;-1]);

%% Centralized solution 

%Combine all matrices into one
Hc = blkdiag(H1,H2,H3,H4,zeros(4,4));
hc = [h1, h2, h3, h4,zeros(1,4)];
A_eq_c1 = [blkdiag(A_eq_1,A_eq_2,A_eq_3,A_eq_4)];
A_eq_c2 = [-eye(4);-eye(4);-eye(4);-eye(4)];
A_eq_c = [A_eq_c1, A_eq_c2];
b_eq_c = [b_eq_1;b_eq_2;b_eq_3;b_eq_4];
A_ineq_c = [blkdiag(A_ineq,A_ineq,A_ineq,A_ineq),zeros(80,4)];
b_ineq_c = ones(2*4*dim.nu*dim.N,1)*umax/Tfinal;

%Solve for optimal solution with quadprog
opts = optimoptions('quadprog', 'Display', 'off');
uopt = quadprog(2*Hc, hc, A_ineq_c, b_ineq_c, A_eq_c,-b_eq_c,[],[],[],opts);

%Final, common state
xfc = uopt(end-3:end);

%% Aircarft state trajectories for centralized case

%Compute the trajectories for each agent
u1 = uopt(1:10);
x1 = [LTI1.x0; predmod1.T*LTI1.x0+predmod1.S*u1];
x11 = x1(1:4:end); x12 = x1(2:4:end); x13 = x1(3:4:end); x14 = x1(4:4:end);
u2 = uopt(11:20);
x2 = [LTI2.x0;predmod2.T*LTI2.x0+predmod2.S*u2];
x21 = x2(1:4:end); x22 = x2(2:4:end); x23 = x2(3:4:end); x24 = x2(4:4:end);
u3 = uopt(21:30);
x3 = [LTI3.x0;predmod3.T*LTI3.x0+predmod3.S*u3];
x31 = x3(1:4:end); x32 = x3(2:4:end); x33 = x3(3:4:end); x34 = x3(4:4:end);
u4 = uopt(31:40);
x4 = [LTI4.x0; predmod4.T*LTI4.x0+predmod4.S*u4];
x41 = x4(1:4:end); x42 = x4(2:4:end); x43 = x4(3:4:end); x44 = x4(4:4:end);
T = [0,1,2,3,4,5];

% Plot the trajectories
figure('Position', [10 50 1200 600])
t = tiledlayout(2,2);
ax1 = nexttile;
plot(T,x11,'Linewidth', 1.2)
hold on
plot(T,x12,'Linewidth', 1.2)
plot(T,x13,'Linewidth', 1.2)
plot(T,x14,'Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax1,'Aircraft 1')
ax1.FontSize = 14;
ax2 = nexttile;
plot(T,x21,'Linewidth', 1.2)
hold on
plot(T,x22,'Linewidth', 1.2)
plot(T,x23,'Linewidth', 1.2)
plot(T,x24,'Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax2,'Aircraft 2')
ax2.FontSize = 14;
ax3 = nexttile;
plot(T,x31,'Linewidth', 1.2)
hold on
plot(T,x32,'Linewidth', 1.2)
plot(T,x33,'Linewidth', 1.2)
plot(T,x34,'Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax3,'Aircraft 3')
ax3.FontSize = 14;
ax4 = nexttile;
plot(T,x41,'Linewidth', 1.2)
hold on
plot(T,x42,'Linewidth', 1.2)
plot(T,x43,'Linewidth', 1.2)
plot(T,x44,'Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax4,'Aircraft 4')
set(gca, 'fontsize', 15)

%% PROBLEM 1 %%
%% 1a. Dual problem 

%Addapt matrices for the dual problem
H1_d = blkdiag(H1,zeros(4,4));
H2_d = blkdiag(H2,zeros(4,4));
H3_d = blkdiag(H3,zeros(4,4));
H4_d = blkdiag(H4,zeros(4,4));
A_eq_1_d = [A_eq_1, -eye(4)];
A_eq_2_d = [A_eq_2, -eye(4)];
A_eq_3_d = [A_eq_3, -eye(4)];
A_eq_4_d = [A_eq_4, -eye(4)];
A_ineq_d =[A_ineq, zeros(20,4)];

%% 1a. Individual optimization problems

%initialize
mu12 = zeros(4,1); 
mu23 = zeros(4,1); 
mu34 = zeros(4,1);
mu41 = zeros(4,1);

%TODO: Choose step size and iterations
a0_ = 0.5;
%a0_ = [0.5,1,2,5,10];

iter = 100;

for r = 1:size(a0_,2)
    mu12 = zeros(4,1); 
    mu23 = zeros(4,1); 
    mu34 = zeros(4,1);
    mu41 = zeros(4,1);
    for i = 1:iter
        a = a0_(r); %a0_(r)/i; %TODO: Choose variable or constant step size
        opts = optimoptions('quadprog', 'Display', 'off');
    %Problem 1
        h1_d = [h1, mu12'-mu41'];
        uopt1 = quadprog(2*H1_d, h1_d, A_ineq_d, b_ineq, A_eq_1_d,-b_eq_1,[],[],[],opts);
    %Problem 2
        h2_d = [h2, mu23'-mu12'];
        uopt2 = quadprog(2*H2_d, h2_d, A_ineq_d, b_ineq, A_eq_2_d,-b_eq_2,[],[],[],opts);
    %Problem 3
        h3_d = [h3, mu34'-mu23'];
        uopt3 = quadprog(2*H3_d, h3_d, A_ineq_d, b_ineq, A_eq_3_d,-b_eq_3,[],[],[],opts);
    %Problem 4
        h4_d = [h4, mu41'-mu34'];
        uopt4 = quadprog(2*H4_d, h4_d, A_ineq_d, b_ineq, A_eq_4_d,-b_eq_4,[],[],[],opts);
    %mu update
        mu12 = mu12 + a*(uopt1(end-3:end)-uopt2(end-3:end));%/norm(uopt1(end-3:end)-uopt2(end-3:end)); %TODO: Activate for normalized gradient
        mu23 = mu23 + a*(uopt2(end-3:end)-uopt3(end-3:end));%/norm(uopt2(end-3:end)-uopt3(end-3:end));
        mu34 = mu34 + a*(uopt3(end-3:end)-uopt4(end-3:end));%/norm(uopt3(end-3:end)-uopt4(end-3:end));
        mu41 = mu41 + a*(uopt4(end-3:end)-uopt1(end-3:end));%/norm(uopt4(end-3:end)-uopt1(end-3:end));
    %Save variables for convergence analysis
        xf_1{r}(:,i) = uopt1(end-3:end);
        xf_2{r}(:,i) = uopt2(end-3:end);
        xf_3{r}(:,i) = uopt3(end-3:end);
        xf_4{r}(:,i) = uopt4(end-3:end);
        mu_save(:,i) = [mu12;mu23;mu34;mu41];
    end
    %Evaluate the error
    for i = 1:iter
    xf1_err{r}(i) = norm(xf_1{r}(:,i)-xfc,2)/norm(xfc,2);
    xf2_err{r}(i) = norm(xf_2{r}(:,i)-xfc,2)/norm(xfc,2);
    xf3_err{r}(i) = norm(xf_3{r}(:,i)-xfc,2)/norm(xfc,2);
    xf4_err{r}(i) = norm(xf_4{r}(:,i)-xfc,2)/norm(xfc,2);
    xf_err_dual{r}(i) = 1/4*(xf1_err{r}(i)+xf2_err{r}(i)+xf3_err{r}(i)+xf4_err{r}(i));
    end
end

%% 1a. Error sequence plot for all agents

figure('Position', [10 50 1200 600])
hold on
for r = 1:size(a0_,2)
    plot(xf1_err{r},'Linewidth', 1.2)
    hold on
    plot(xf2_err{r},'Linewidth', 1.2)
    plot(xf3_err{r},'Linewidth', 1.2)
    plot(xf4_err{r},'Linewidth', 1.2)
end
xlabel('Iterations k')
ylabel('$\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}$','Interpreter','latex')
legend('Aircraft 1', 'Aircraft 2','Aircraft 3','Aircraft 4')
set(gca, 'fontsize', 15)

%% 1b. Error sequence plot for different stepsizes (mean of agents)

figure('Position', [10 50 1200 600])
hold on
for r = 1:size(a0_,2)
    plot(xf_err_dual{r},'Linewidth', 1.2)
end
xlabel('Iterations k')
ylabel('$\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}$','Interpreter','latex')
legend('\alpha = 0.5', '\alpha = 1','\alpha = 2','\alpha = 5','\alpha = 10')
set(gca, 'fontsize', 15)

%% 1a. Plot the trajectories in comparison

x1d = [LTI1.x0; predmod1.T*LTI1.x0+predmod1.S*uopt1(1:end-4)];
x11d = x1d(1:4:end); x12d = x1d(2:4:end); x13d = x1d(3:4:end); x14d = x1d(4:4:end);
x2d = [LTI2.x0;predmod2.T*LTI2.x0+predmod2.S*uopt2(1:end-4)];
x21d = x2d(1:4:end); x22d = x2d(2:4:end); x23d = x2d(3:4:end); x24d = x2d(4:4:end);
x3d = [LTI3.x0;predmod3.T*LTI3.x0+predmod3.S*uopt3(1:end-4)];
x31d = x3d(1:4:end); x32d = x3d(2:4:end); x33d = x3d(3:4:end); x34d = x3d(4:4:end);
x4d = [LTI4.x0; predmod4.T*LTI4.x0+predmod4.S*uopt4(1:end-4)];
x41d = x4(1:4:end); x42d = x4(2:4:end); x43d = x4(3:4:end); x44d = x4(4:4:end);
T = [0,1,2,3,4,5];

figure('Position', [10 50 1200 600])
t = tiledlayout(2,2);
ax1 = nexttile;
plot(T,x11,'Linewidth', 1.2)
hold on
plot(T,x12,'Linewidth', 1.2)
plot(T,x13,'Linewidth', 1.2)
plot(T,x14,'Linewidth', 1.2)
plot(T,x11d,'c--','Linewidth', 1.2)
plot(T,x12d,'c--','Linewidth', 1.2)
plot(T,x13d,'c--','Linewidth', 1.2)
plot(T,x14d,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax1,'Aircraft 1')
ax1.FontSize = 14;
ax2 = nexttile;
plot(T,x21,'Linewidth', 1.2)
hold on
plot(T,x22,'Linewidth', 1.2)
plot(T,x23,'Linewidth', 1.2)
plot(T,x24,'Linewidth', 1.2)
plot(T,x21d,'c--','Linewidth', 1.2)
plot(T,x22d,'c--','Linewidth', 1.2)
plot(T,x23d,'c--','Linewidth', 1.2)
plot(T,x24d,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax2,'Aircraft 2')
ax2.FontSize = 14;
ax3 = nexttile;
plot(T,x31,'Linewidth', 1.2)
hold on
plot(T,x32,'Linewidth', 1.2)
plot(T,x33,'Linewidth', 1.2)
plot(T,x34,'Linewidth', 1.2)
plot(T,x31d,'c--','Linewidth', 1.2)
plot(T,x32d,'c--','Linewidth', 1.2)
plot(T,x33d,'c--','Linewidth', 1.2)
plot(T,x34d,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax3,'Aircraft 3')
ax3.FontSize = 14;
ax4 = nexttile;
plot(T,x41,'Linewidth', 1.2)
hold on
plot(T,x42,'Linewidth', 1.2)
plot(T,x43,'Linewidth', 1.2)
plot(T,x44,'Linewidth', 1.2)
plot(T,x41d,'c--','Linewidth', 1.2)
plot(T,x42d,'c--','Linewidth', 1.2)
plot(T,x43d,'c--','Linewidth', 1.2)
plot(T,x44d,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax4,'Aircraft 4')
ax4.FontSize = 14;

%% 1c. Accelerated gradient 

%A-Matrix for Nesterovs accelerated gradient
A_eq_acc = [A_eq_1, -A_eq_2, zeros(size(A_eq_1)), zeros(size(A_eq_1)); ...
    zeros(size(A_eq_1)), A_eq_2, - A_eq_3, zeros(size(A_eq_1)); ...
    zeros(size(A_eq_1)),zeros(size(A_eq_1)), A_eq_3, - A_eq_4;...
    -A_eq_1, zeros(size(A_eq_1)), zeros(size(A_eq_1)), A_eq_4];
A_ineq_acc = blkdiag(A_ineq, A_ineq, A_ineq, A_ineq);
A_acc = [A_eq_acc; A_ineq_acc];

%b-Matrix for Nesterovs accelerated gradient
b_eq_acc = [b_eq_2 - b_eq_1; b_eq_3 - b_eq_2; b_eq_4 - b_eq_3; b_eq_1 - b_eq_4];
b_ineq_acc = [b_ineq;b_ineq;b_ineq;b_ineq];
b_acc = [b_eq_acc; b_ineq_acc];

%combined H-marix
H_acc = blkdiag(H1,H2,H3,H4);
h_acc = [h1, h2, h3, h4];

%Lipschitz constant
L = norm(0.5*A_acc*inv(H_acc)*A_acc',2);

%Initialize lambdas and mus
l12 = zeros(4,1); l23 = zeros(4,1); l34 = zeros(4,1); l41 = zeros(4,1);
mu1 = zeros(20,1); mu2 = zeros(20,1); mu3 = zeros(20,1); mu4 = zeros(20,1);

% Initialize z's
z0 = [l12; l23; l34; l41; mu1; mu2; mu3; mu4];
z1 = [l12; l23; l34; l41; mu1; mu2; mu3; mu4];
z(:,1) = z0;
z(:,2) = z1;

for i = 1:100
    opts = optimoptions('quadprog', 'Display', 'off');
%Problem 1
    h1_acc = (h1 +(l12'-l41')*A_eq_1 + mu1'*A_ineq);
    uopt1 = quadprog(2*H1, h1_acc,[],[],[],[],[],[],[],opts);
    xf1(:,i) = A_eq_1*uopt1 + b_eq_1;
%Problem 2
    h2_acc = (h2 +(l23'-l12')*A_eq_2 + mu2'*A_ineq);
    uopt2 = quadprog(2*H2, h2_acc,[],[],[],[],[],[],[],opts);
    xf2(:,i) = A_eq_2*uopt2 + b_eq_2;
%Problem 3
    h3_acc = (h3 +(l34'-l23')*A_eq_3 + mu3'*A_ineq);
    uopt3 = quadprog(2*H3, h3_acc,[],[],[],[],[],[],[],opts);
    xf3(:,i) = A_eq_3*uopt3 + b_eq_3;
%Problem 4
    h4_acc = (h4 +(l41'-l34')*A_eq_4 + mu4'*A_ineq);
    uopt4 = quadprog(2*H4, h4_acc,[],[],[],[],[],[],[],opts);
    xf4(:,i) = A_eq_4*uopt4 + b_eq_4;

%Update
    vk = z(:,i+1) + .75*(i-1)/(i+2)*(z(:,i+1)-z(:,i));
    z(:,i+2) = vk - 1/L*(1/2*A_acc/H_acc*(A_acc'*vk+h_acc') + b_acc);
    ls(:,i+2) = z(1:16,i+2); % ls = [l-1 , l0, l1, l2, ...]
    mus(:,i+2) = max(z(17:end,i+2), 0); % mus = [mu-1, mu0, mu1, mu2, ...]
    z(:,i+2) = [ls(:,i+2); mus(:,i+2)];

    l12 = ls(1:4,i+2); l23 = ls(5:8,i+2); l34 = ls(9:12,i+2); l41 = ls(13:16,i+2);   
    mu1 = mus(1:20,i+2); mu2 = mus(21:40,i+2); mu3 = mus(41:60,i+2); mu4 = mus(61:80,i+2);
end

%% 1c. Error sequence plot compared with 1a.

for i = 1:100
    xf1_err_acc(i) = norm(xf1(:,i)-xfc,2)/norm(xfc,2);
    xf2_err_acc(i) = norm(xf2(:,i)-xfc,2)/norm(xfc,2);
    xf3_err_acc(i) = norm(xf3(:,i)-xfc,2)/norm(xfc,2);
    xf4_err_acc(i) = norm(xf4(:,i)-xfc,2)/norm(xfc,2);
end

figure('Position', [10 50 1200 600])
plot(xf1_err{1},'--','Linewidth', 1.2)
hold on
plot(xf2_err{1},'--','Linewidth', 1.2)
plot(xf3_err{1},'--','Linewidth', 1.2)
plot(xf4_err{1},'--','Linewidth', 1.2)
plot(xf1_err_acc,'Linewidth', 1.2)
plot(xf2_err_acc,'Linewidth', 1.2)
plot(xf3_err_acc,'Linewidth', 1.2)
plot(xf4_err_acc,'Linewidth', 1.2)
%plot(xf_err_admm,'r','Linewidth', 1.5)
xlabel('Iterations k')
ylabel('$\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}$','Interpreter','latex')
legend('Aircraft 1 Dual', 'Aircraft 2 Dual','Aircraft 3 Dual','Aircraft 4 Dual','Aircraft 1 Accelerated', 'Aircraft 2 Accelerated','Aircraft 3 Accelerated','Aircraft 4 Accelerated')%, 'ADMM')
set(gca, 'fontsize', 15)
%set(gca,'ColorOrder',Customcolor)

%% 1c. Plot the trajectories in comparison

x1a = [LTI1.x0; predmod1.T*LTI1.x0+predmod1.S*uopt1];
x11a = x1a(1:4:end); x12a = x1a(2:4:end); x13a = x1a(3:4:end); x14a = x1a(4:4:end);
x2a = [LTI2.x0;predmod2.T*LTI2.x0+predmod2.S*uopt2];
x21a = x2a(1:4:end); x22a = x2a(2:4:end); x23a = x2a(3:4:end); x24a = x2a(4:4:end);
x3a = [LTI3.x0;predmod3.T*LTI3.x0+predmod3.S*uopt3];
x31a = x3a(1:4:end); x32a = x3a(2:4:end); x33a = x3a(3:4:end); x34a = x3a(4:4:end);
x4a = [LTI4.x0; predmod4.T*LTI4.x0+predmod4.S*uopt4];
x41a = x4a(1:4:end); x42a = x4a(2:4:end); x43a = x4a(3:4:end); x44a = x4a(4:4:end);
T = [0,1,2,3,4,5];

figure('Position', [10 50 1200 600])
t = tiledlayout(2,2);
ax1 = nexttile;
plot(T,x11,'Linewidth', 1.2)
hold on
plot(T,x12,'Linewidth', 1.2)
plot(T,x13,'Linewidth', 1.2)
plot(T,x14,'Linewidth', 1.2)
plot(T,x11a,'c--','Linewidth', 1.2)
plot(T,x12a,'c--','Linewidth', 1.2)
plot(T,x13a,'c--','Linewidth', 1.2)
plot(T,x14a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax1,'Aircraft 1')
ax1.FontSize = 14;
ax2 = nexttile;
plot(T,x21,'Linewidth', 1.2)
hold on
plot(T,x22,'Linewidth', 1.2)
plot(T,x23,'Linewidth', 1.2)
plot(T,x24,'Linewidth', 1.2)
plot(T,x21a,'c--','Linewidth', 1.2)
plot(T,x22a,'c--','Linewidth', 1.2)
plot(T,x23a,'c--','Linewidth', 1.2)
plot(T,x24a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax2,'Aircraft 2')
ax2.FontSize = 14;
ax3 = nexttile;
plot(T,x31,'Linewidth', 1.2)
hold on
plot(T,x32,'Linewidth', 1.2)
plot(T,x33,'Linewidth', 1.2)
plot(T,x34,'Linewidth', 1.2)
plot(T,x31a,'c--','Linewidth', 1.2)
plot(T,x32a,'c--','Linewidth', 1.2)
plot(T,x33a,'c--','Linewidth', 1.2)
plot(T,x34a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax3,'Aircraft 3')
ax3.FontSize = 14;
ax4 = nexttile;
plot(T,x41,'Linewidth', 1.2)
hold on
plot(T,x42,'Linewidth', 1.2)
plot(T,x43,'Linewidth', 1.2)
plot(T,x44,'Linewidth', 1.2)
plot(T,x41a,'c--','Linewidth', 1.2)
plot(T,x42a,'c--','Linewidth', 1.2)
plot(T,x43a,'c--','Linewidth', 1.2)
plot(T,x44a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax4,'Aircraft 4')
ax4.FontSize = 14;
%% 1d. Consensus approach

%initialize
W = [0.75, 0.25, 0.00, 0.00; ...
     0.25, 0.50, 0.25, 0.00; ...
     0.00, 0.25, 0.50, 0.25; ...
     0.00, 0.00, 0.25, 0.75];
 
%TODO: Chose sequence of phi's to analyze and iterations
phi_ = [0,1,2,5,10,20,50];

iter = 50;

for r = 1:size(phi_,2)
    phi = phi_(r)
    %Compute consensus matrix after phi updates
    Wphi = W ^ phi;
    
    %Initialize
    clear xf
    xf{1} = zeros(4,iter+1);
    xf{2} = zeros(4,iter+1);
    xf{3} = zeros(4,iter+1);
    xf{4} = zeros(4,iter+1); 

    for k = 1:iter
        a = 0.5/(2*(k+1)); %TODO: Chose variable or constant step size
        opts = optimoptions('quadprog', 'Display', 'off');
    %Problem 1
        [uopt1,~,~,~,lambda] = quadprog(2*H1, h1, A_ineq, b_ineq, A_eq_1, xf{1}(:,k)-b_eq_1,[],[],[],opts);
        g{1}(:,k) = lambda.eqlin;
    %Problem 2
        [uopt2,~,~,~,lambda] = quadprog(2*H2, h2, A_ineq, b_ineq, A_eq_2,xf{2}(:,k)-b_eq_2,[],[],[],opts);
        g{2}(:,k) = lambda.eqlin;
    %Problem 3
        [uopt3,~,~,~,lambda] = quadprog(2*H3, h3, A_ineq, b_ineq, A_eq_3,xf{3}(:,k)-b_eq_3,[],[],[],opts);
        g{3}(:,k) = lambda.eqlin;
    %Problem 4
        [uopt4,~,~,~,lambda] = quadprog(2*H4, h4, A_ineq, b_ineq, A_eq_4,xf{4}(:,k)-b_eq_4,[],[],[],opts);
        g{4}(:,k) = lambda.eqlin;
      
    %consensus update 
        for i = 1:4
            for j = 1:4
                xf{i}(:, k + 1) = xf{i}(:, k + 1) + Wphi(i, j) * (xf{j}(:, k) + a.* g{j}(:, k));
            end
        end
    end
    
    %Error sequence
    for i = 1:iter
        xf1_err_c(i) = norm(xf{1}(:,i)-xfc,2)/norm(xfc,2);
        xf2_err_c(i) = norm(xf{2}(:,i)-xfc,2)/norm(xfc,2);
        xf3_err_c(i) = norm(xf{3}(:,i)-xfc,2)/norm(xfc,2);
        xf4_err_c(i) = norm(xf{4}(:,i)-xfc,2)/norm(xfc,2);
        xf_err_c{r}(i) = (xf1_err_c(i) + xf1_err_c(i) + xf1_err_c(i) + xf1_err_c(i))/4;
    end
end 
%% Plot error sequence for different consensus parameters

it = 1:1:iter;

figure('Position', [10 50 1200 600])
hold on
for i = 1:7
    plot(it,log10(xf_err_c{i}),'Linewidth', 1.2)
end
xlabel('Iterations k')
ylabel('$log_{10}\left(\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}\right)$','Interpreter','latex')
legend('\Phi = 0','\Phi = 1','\Phi = 2','\Phi = 5','\Phi = 10', '\Phi = 20','\Phi = 50')
set(gca, 'fontsize', 15)

%% PROBLEM 2 %%
%% 2a. ADMM

%initialize
xf0 = [0;0;0;0];
y0 = [0;0;0;0];
rho = 2;
iter = 100;

xf = zeros(N,iter);
xf(:,1) = xf0;
y1 = zeros(dim.nx, iter);
y1(:,1) = y0;
y2 = zeros(dim.nx, iter);
y2(:,1) = y0;
y3 = zeros(dim.nx, iter);
y3(:,1) = y0;
y4 = zeros(dim.nx, iter);
y4(:,1) = y0;
uopt1_ = zeros(iter, dim.nu*dim.N);
uopt2_ = zeros(iter, dim.nu*dim.N);
uopt3_ = zeros(iter, dim.nu*dim.N);
uopt4_ = zeros(iter, dim.nu*dim.N);
ubar = zeros(N,iter);
ubar(:,1) = y0;

%Solve algorithm with quadprog
rho_ = [0.1,0.5,1,2,5,10,20];
rho_ = 10;
for r = 1:size(rho_,2)
    rho = rho_(r)
    clear xf
    xf{r}(:,1) = xf0;
    for k = 1:iter
        %1) Update x's
        %u1(k+1)
        uopt1 = quadprog(2*(H1 + rho/2.*A_eq_1'*A_eq_1),(h1+y1(:,k)'*A_eq_1+rho*(b_eq_1 - xf{r}(:,k))'*A_eq_1),A_ineq, b_ineq,[],[],[],[],[],opts);
        uopt1_(k,:) = uopt1;
        %u2(k+1)
        uopt2 = quadprog(2*(H2 + rho/2.*A_eq_2'*A_eq_2),(h2+y2(:,k)'*A_eq_2+rho*(b_eq_2 - xf{r}(:,k))'*A_eq_2),A_ineq, b_ineq,[],[],[],[],[],opts);
        uopt2_(k,:) = uopt2;
        %u3(k+1)
        uopt3 = quadprog(2*(H3 + rho/2.*A_eq_3'*A_eq_3),(h3+y3(:,k)'*A_eq_3+rho*(b_eq_3 - xf{r}(:,k))'*A_eq_3),A_ineq, b_ineq,[],[],[],[],[],opts);
        uopt3_(k,:) = uopt3;
        %u4(k+1)
        uopt4 = quadprog(2*(H4 + rho/2.*A_eq_4'*A_eq_4),(h4+y4(:,k)'*A_eq_4+rho*(b_eq_4 - xf{r}(:,k))'*A_eq_4),A_ineq, b_ineq,[],[],[],[],[],opts);
        uopt4_(k,:) = uopt4;
        %2) Update xf
        xf{r}(:,k+1) = 1/N*((A_eq_1*uopt1 + b_eq_1 + 1/rho*y1(:,k))+(A_eq_2*uopt2 + b_eq_2 + 1/rho*y2(:,k))+(A_eq_3*uopt3 + b_eq_3 + 1/rho*y3(:,k))+(A_eq_4*uopt4 + b_eq_4 + 1/rho*y4(:,k)));
        %3) Update y's
        y1(:,k+1) = y1(:,k) + rho*(A_eq_1*uopt1 + b_eq_1 - xf{r}(:,k+1));
        y2(:,k+1) = y2(:,k) + rho*(A_eq_2*uopt2 + b_eq_2 - xf{r}(:,k+1));
        y3(:,k+1) = y3(:,k) + rho*(A_eq_3*uopt3 + b_eq_3 - xf{r}(:,k+1));
        y4(:,k+1) = y4(:,k) + rho*(A_eq_4*uopt4 + b_eq_4 - xf{r}(:,k+1));
    end
    for i = 1:100
        xf_err_admm{r}(i) = norm(xf{r}(:,i)-xfc,2)/norm(xfc,2);
    end
end 
%% 2a. Error sequence plot (also for 2b. for muliple rho's)
figure('Position', [10 50 1200 600])
hold on
for r = 1:size(rho_,2)
    plot(xf_err_admm{r},'Linewidth', 1.2)
end
xlabel('Iterations k')
ylabel('$\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}$','Interpreter','latex')
legend('\rho = 0.1','\rho = 0.5','\rho = 1','\rho = 2','\rho = 5','\rho = 10','\rho = 20')
set(gca, 'fontsize', 15)

%% 2a. Compare convergnece to other methods (accelerated gradient and dual approach)

xf_err_acc = 1/4*(xf1_err_acc+xf2_err_acc+xf3_err_acc+xf4_err_acc);
figure('Position', [10 50 1200 600])
plot(xf_err_dual{1},'Linewidth', 1.2)
hold on
plot(xf_err_acc,'Linewidth', 1.2)
plot(xf_err_admm{1},'Linewidth', 1.2)
xlabel('Iterations k')
ylabel('$\frac{||x_f^i(k)-x_f^*||}{||x_f^*||}$','Interpreter','latex')
legend('Dual', 'Accelerated Gradient','ADMM')
set(gca, 'fontsize', 15)

%% 2a. Plot the trajectories in comparison

x1a = [LTI1.x0; predmod1.T*LTI1.x0+predmod1.S*uopt1];
x11a = x1a(1:4:end); x12a = x1a(2:4:end); x13a = x1a(3:4:end); x14a = x1a(4:4:end);
x2a = [LTI2.x0;predmod2.T*LTI2.x0+predmod2.S*uopt2];
x21a = x2a(1:4:end); x22a = x2a(2:4:end); x23a = x2a(3:4:end); x24a = x2a(4:4:end);
x3a = [LTI3.x0;predmod3.T*LTI3.x0+predmod3.S*uopt3];
x31a = x3a(1:4:end); x32a = x3a(2:4:end); x33a = x3a(3:4:end); x34a = x3a(4:4:end);
x4a = [LTI4.x0; predmod4.T*LTI4.x0+predmod4.S*uopt4];
x41a = x4a(1:4:end); x42a = x4a(2:4:end); x43a = x4a(3:4:end); x44a = x4a(4:4:end);
T = [0,1,2,3,4,5];

figure('Position', [10 50 1200 600])
t = tiledlayout(2,2);
ax1 = nexttile;
plot(T,x11,'Linewidth', 1.2)
hold on
plot(T,x12,'Linewidth', 1.2)
plot(T,x13,'Linewidth', 1.2)
plot(T,x14,'Linewidth', 1.2)
plot(T,x11a,'c--','Linewidth', 1.2)
plot(T,x12a,'c--','Linewidth', 1.2)
plot(T,x13a,'c--','Linewidth', 1.2)
plot(T,x14a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax1,'Aircraft 1')
ax1.FontSize = 14;
ax2 = nexttile;
plot(T,x21,'Linewidth', 1.2)
hold on
plot(T,x22,'Linewidth', 1.2)
plot(T,x23,'Linewidth', 1.2)
plot(T,x24,'Linewidth', 1.2)
plot(T,x21a,'c--','Linewidth', 1.2)
plot(T,x22a,'c--','Linewidth', 1.2)
plot(T,x23a,'c--','Linewidth', 1.2)
plot(T,x24a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax2,'Aircraft 2')
ax2.FontSize = 14;
ax3 = nexttile;
plot(T,x31,'Linewidth', 1.2)
hold on
plot(T,x32,'Linewidth', 1.2)
plot(T,x33,'Linewidth', 1.2)
plot(T,x34,'Linewidth', 1.2)
plot(T,x31a,'c--','Linewidth', 1.2)
plot(T,x32a,'c--','Linewidth', 1.2)
plot(T,x33a,'c--','Linewidth', 1.2)
plot(T,x34a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax3,'Aircraft 3')
ax3.FontSize = 14;
ax4 = nexttile;
plot(T,x41,'Linewidth', 1.2)
hold on
plot(T,x42,'Linewidth', 1.2)
plot(T,x43,'Linewidth', 1.2)
plot(T,x44,'Linewidth', 1.2)
plot(T,x41a,'c--','Linewidth', 1.2)
plot(T,x42a,'c--','Linewidth', 1.2)
plot(T,x43a,'c--','Linewidth', 1.2)
plot(T,x44a,'c--','Linewidth', 1.2)
legend('x_1','x_2','x_3','x_4',NumColumns = 2)
title(ax4,'Aircraft 4')
ax4.FontSize = 14;