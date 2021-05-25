%% Initially given
%Student number: 4677498
clear all
clc
a = 4;
b = 7;
c = 8;
A = [a b+0.5; 0 -c];
B = [0;1];
% Compute controller through pole placement
p = [-2 -3];
K = place(A,B,p);

%% Question 3
%% 3.1 Maximum Number of packet dropouts
syms h hl
Ai = expm(A*hl)-expm(A*(hl-h))*(expm(A*h)-eye(2))*inv(A)*B*K; %Symbolic expression for the switched system matrix
Q = eye(2);
h_eval = linspace(0.01,0.3,31)
for i = 1:size(h_eval,2)
    % System matrices evaluated for increasing number of package dropouts
    h = h_eval(i);
    hl = h;
    A01 = eval(Ai);
    hl = 2*h;
    A11 = eval(Ai);
    hl = 3*h;
    A21 = eval(Ai);
    hl = 4*h;
    A31 = eval(Ai);
    hl = 5*h;
    A41 = eval(Ai);
    %Solve set of LMIs, adding more switched system matrices if solution is
    %found
    cvx_begin sdp 
    variable P(2,2)
    subject to 
        A01'*P*A01 - P <= -Q;
        P >= eye(2);
    cvx_end
    c = -1;
    if (strcmp(cvx_status, 'Solved'))
        c = 0;
        cvx_begin sdp 
        variable P(2,2)
        subject to 
            A01'*P*A01 - P <= -Q;
            A11'*P*A11 - P <= -Q;
            P >= eye(2);
        cvx_end
        if (strcmp(cvx_status, 'Solved'))
            c = 1;
            cvx_begin sdp 
            variable P(2,2)
            subject to 
                A01'*P*A01 - P <= -Q;
                A11'*P*A11 - P <= -Q;
                A21'*P*A21 - P <= -Q;
            P >= eye(2);
            cvx_end
            if (strcmp(cvx_status, 'Solved'))
                c = 2;
                cvx_begin sdp 
                variable P(2,2)
                subject to 
                    A01'*P*A01 - P <= -Q;
                    A11'*P*A11 - P <= -Q;
                    A21'*P*A21 - P <= -Q;
                    A31'*P*A31 - P <= -Q;
                    P >= eye(2);
                cvx_end
                if (strcmp(cvx_status, 'Solved'))
                    c= 3;
                    cvx_begin sdp 
                    variable P(2,2)
                    subject to 
                        A01'*P*A01 - P <= -Q;
                        A11'*P*A11 - P <= -Q;
                        A21'*P*A21 - P <= -Q;
                        A31'*P*A31 - P <= -Q;
                        A41'*P*A41 - P <= -Q;
                    P >= eye(2);
                    cvx_end
                    if (strcmp(cvx_status, 'Solved'))
                        c = 4;
                    end
                end
            end
        end
    end
    %Store count
max_drop(i) = c; 
end

%% Plot the results
figure(1)
plot(h_eval(1:end),max_drop(1:end),'ro-')
xlabel('Sampling time h')
ylabel('\delta')
title('Maximum Number of Packet Dropouts Under Which the NCS Remains Stable')

%% 3.2 Upper bound on p for Bernoulli stability

h = 0.1;
F = expm(A*h);
G = (expm(A*h)-eye(2))/A*B;

A0 = F-G*K;
A1 = F;
Q = eye(2);

%Solve LMI for incerasing probability, stop if no longer solvable
for p = linspace(0,0.5,501)
    cvx_begin sdp 
        variable P(2,2)
        subject to 
            P - (1-p)*A0'*P*A0 - p*A1'*P*A1 >= Q;
            P >= eye(2);
    cvx_end
    if strcmp(cvx_status, 'Infeasible')
        break
    end
end

% Problem stable until p= 0.103

%% 3.3 Verify results by numerical simulation 

%Simulate the system for different inital conditions and probabilities
N = 100;
x0 = [1;1];
xk = zeros(2,N+1);
xk(:,1) = x0;
p_star = 0.3; %Change probability 

%Simulate a few times and plot the results
for i=1:N
    c = rand;
    if c < p_star
        xk(:,i+1) = A1*xk(:,i);
    else
        xk(:,i+1) = A0*xk(:,i);
    end
end

figure(2)
hold on 

title('Trajectory for multiple simulations with probabiliy p < p*')
subplot(2,2,1); 
plot(xk(1,:),xk(2,:))
xlabel('x1')
ylabel('x2')

xk = zeros(2,N+1);
xk(:,1) = x0;
for i=1:N
    c = rand;
    if c < p_star
        xk(:,i+1) = A1*xk(:,i);
    else
        xk(:,i+1) = A0*xk(:,i);
    end
end

subplot(2,2,3); 
plot(xk(1,:),xk(2,:))
xlabel('x1')
ylabel('x2')

xk = zeros(2,N+1);
xk(:,1) = x0;
for i=1:N
    c = rand;
    if c < p_star
        xk(:,i+1) = A1*xk(:,i);
    else
        xk(:,i+1) = A0*xk(:,i);
    end
end

subplot(2,2,2); 
plot(xk(1,:),xk(2,:))
xlabel('x1')
ylabel('x2')

xk = zeros(2,N+1);
xk(:,1) = x0;
for i=1:N
    c = rand;
    if c < p_star
        xk(:,i+1) = A1*xk(:,i);
    else
        xk(:,i+1) = A0*xk(:,i);
    end
end

subplot(2,2,4); 
plot(xk(1,:),xk(2,:))
xlabel('x1')
ylabel('x2')

%% 3.3 Look at expected value 
N = 10000;
x0 = [-100;-100]; %Test for different initial conditions
xk = zeros(2,N+1);
xk(:,1) = x0;
p_star = 0.1; %Change probability 
%Simulate system trajectory 1000 times for k--> inf
for k = 1:1000;
    for i=1:N
        c = rand;
        if c < p_star
            xk(:,i+1) = A1*xk(:,i);
        else
            xk(:,i+1) = A0*xk(:,i);
        end
    end
    x(:,k) = xk(:,end);
end
%Compute the mean
x_mean_1 = mean(x(1,:));
x_mean_2 = mean(x(2,:));

%% 3.4 Upper Bound for Bernoully stability with ASS theory

h = 0.1;
F = expm(A*h);
G = (expm(A*h)-eye(2))/A*B;
A0 = F-G*K;
A1 = F;
p_eval = linspace(0.01,0.2,20);

for i = 1:size(p_eval,2);
    p = p_eval(i);
    p_stab(i) = 0;
    for lambda = linspace(0.1,1-eps,10); %Vary lambda between 0 and 1-eps
        L = (1/lambda^(1-p))^(1/p) - eps; % Line search problem
        cvx_begin sdp 
            variable P(2,2)
            subject to 
                A0'*P*A0 <= lambda*P
                A1'*P*A1 <= L*P
                P >= eye(2);
        cvx_end
        if strcmp(cvx_status, 'Solved') %Stop if one solution could be found
            p_stab(i) = 1; %Change array entry to 1 to indicate that the system is stable
            break
        end
    end
end

%% 3.5 Stability for Gilbert-Elliot Process

h = 0.1;
F = expm(A*h);
G = (expm(A*h)-eye(2))/A*B;
A0 = F-G*K;
A1 = F;

%Evaluate stability for different combinations of p and q
[peval, qeval] = meshgrid(linspace(0, 1, 21), linspace(0, 1, 21));
stab = zeros(3,1);
unstab = zeros(3,1);
for i = 1:21*21
    p = peval(i);
    q = qeval(i);
    cvx_begin sdp 
        variable P0(2,2)
        variable P1(2,2)
        subject to 
            P0 - p*A0'*P0*A0 - (1-p)*A1'*P1*A1 >= eye(2);
            P1 - q*A1'*P1*A1 - (1-q)*A0'*P0*A0 >= eye(2);
            P0 >= eye(2);
            P1 >= eye(2);
    cvx_end
    if strcmp(cvx_status, 'Solved')
        stab = [stab, [p;q;1]]; %Store stable points
    end
    if strcmp(cvx_status, 'Infeasible')
        unstab = [unstab, [p;q;0]]; %Otherwise store unstable points
    end
end

%% Plot the stable and unstable points
figure(2)
scatter(unstab(1,2:end),unstab(2,2:end),'r')
hold on
scatter(stab(1,2:end),stab(2,2:end),'g')
xlabel('p')
ylabel('q')
title('Stability for combinations of p and q')

%% Question 4
%% 4.1 Jordan form approach 
[Q,J] = jordan(A);
Q = inv(Q);
syms h tau

%System matrices and coefficients from Jordan form Approach
F = expm(A*h);
F0_h = [F, (Q\[0 0; 0 1/4*exp(4*h)]*Q + Q\[-1/8*exp(-8*h) 0; 0 0]*Q)*B; zeros(1,3)];
F1 = [zeros(2,2), Q\[1 0; 0 0]*Q*B; zeros(1,3)];
F2 = [zeros(2,2), Q\[0 0; 0 -1]*Q*B; zeros(1,3)];
G0 = [(Q\[0 0; 0 -1/4]*Q + Q\[1/8 0; 0 0]*Q)*B; 1];
G1 = [Q\[-1 0; 0 0]*Q*B; 0];
G2 = [Q\[0 0; 0 1]*Q*B; 0];
a1 = 1/8*exp(-8*h + 8*tau);
a2 = 1/4*exp(4*h - 4*tau);

%% 4.2 Stability Analysis using LMI's
h_eval = linspace(0, 0.3, 31);
c = 0;
d = 0;
Q4 = eye(3);
K4 = [K, 0]; %Extended controller
clear Abar
clear stab_points2
clear points 
for i = 1:31
    h = h_eval(i);
    for tau = linspace(0,h,31)
        taumax = tau;
        taumin = 0;
        F0 = eval(F0_h);
        %Evaluate at taumax & taumin
        a1max = 1/8*exp(-8*h + 8*taumax);
        a1min = 1/8*exp(-8*h + 8*taumin);
        a2max = 1/4*exp(4*h - 4*taumax);
        a2min = 1/4*exp(4*h - 4*taumin);
        %Possible combinations
        HF1 = F0 + a1max*F1 + a2max*F2;
        HG1 = G0 + a1max*G1 + a2max*G2;
        HF2 = F0 + a1max*F1 + a2min*F2;
        HG2 = G0 + a1max*G1 + a2min*G2;
        HF3 = F0 + a1min*F1 + a2max*F2;
        HG3 = G0 + a1min*G1 + a2max*G2;
        HF4 = F0 + a1min*F1 + a2min*F2;
        HG4 = G0 + a1min*G1 + a2min*G2;
        d = d+1;
        points(:,d) = [h,tau];
        %Solve LMI
        cvx_begin sdp 
            variable P(3,3)
            subject to 
            (HF1 - HG1*K4)'*P*(HF1 - HG1*K4) - P <= -Q4;
            (HF2 - HG2*K4)'*P*(HF2 - HG2*K4) - P <= -Q4;
            (HF3 - HG3*K4)'*P*(HF3 - HG3*K4) - P <= -Q4;
            (HF4 - HG4*K4)'*P*(HF4 - HG4*K4) - P <= -Q4;
            P >= eye(3);
        cvx_end
        %Indicate stability if LMI solved
        if strcmp(cvx_status, 'Solved');
            c = c+1;
            stab_points2(:,c) = [h,tau];
        end   
    end
end
%% Plot the stable and unstable points
figure(3)
hold on
for i = 1:size(points,2)
    plot(points(1,i),points(2,i),'r*');
end
for k = 1:size(stab_points2,2)
    plot(stab_points2(1,k),stab_points2(2,k),'g*');
end
xlabel('Sampling time h')
ylabel('Delay \tau')
title('Stability analysis for different combinations of h and \tau')
%% 4.3 Add more vertices to the system
%Analyze where the vertices have to lie by looking at evolution of a1 and
%a2 
h = 0.1;
tau_eval = linspace(0, 0.1, 11);
for i = 1:11
    taumax = tau_eval(i);
    a1(i) = 1/8*exp(-8*h + 8*taumax);
    a2(i) = 1/4*exp(4*h - 4*taumax);
end

plot(a1,a2);
xlabel('a_1(\tau)')
ylabel('a_2(\tau)')
title('Evolution of a_1 and a_2 as a function of \tau')
%% 4.3 Compute stable and unstable regions
%Same as before now with 6 possible combinations
h_eval = linspace(0, 0.3, 31);
c = 0;
d = 0;
Q4 = eye(3);
K4 = [K, 0];
clear Abar
clear stab_points3
clear points 
for i = 1:31
    h = h_eval(i);
    for tau = linspace(0,h,31)
        taumax = tau;
        taumin = 0;
        F0 = eval(F0_h);
        a1max = 1/8*exp(-8*h + 8*taumax);
        a1med = 1/8*exp(-8*h + 8*taumax/2);
        a1min = 1/8*exp(-8*h + 8*taumin);
        a2max = 1/4*exp(4*h - 4*taumax);
        a2med = 1/4*exp(4*h - 4*taumax/2);
        a2min = 1/4*exp(4*h - 4*taumin);
        HF1 = F0 + a1max*F1 + a2min*F2;
        HG1 = G0 + a1max*G1 + a2min*G2;
        HF2 = F0 + a1min*F1 + a2max*F2;
        HG2 = G0 + a1min*G1 + a2max*G2;
        HF3 = F0 + a1max*F1 + a2med*F2;
        HG3 = G0 + a1max*G1 + a2med*G2;
        HF4 = F0 + a1med*F1 + a2max*F2;
        HG4 = G0 + a1med*G1 + a2max*G2;
        HF5 = F0 + a1min*F1 + a2med*F2;
        HG5 = G0 + a1min*G1 + a2med*G2;
        HF6 = F0 + a1med*F1 + a2min*F2;
        HG6 = G0 + a1med*G1 + a2min*G2;
        d = d+1;
        points(:,d) = [h,tau];
        cvx_begin sdp 
            variable P(3,3)
            subject to 
            (HF1 - HG1*K4)'*P*(HF1 - HG1*K4) - P <= -Q4;
            (HF2 - HG2*K4)'*P*(HF2 - HG2*K4) - P <= -Q4;
            (HF3 - HG3*K4)'*P*(HF3 - HG3*K4) - P <= -Q4;
            (HF4 - HG4*K4)'*P*(HF4 - HG4*K4) - P <= -Q4;
            (HF5 - HG5*K4)'*P*(HF5 - HG5*K4) - P <= -Q4;
            (HF6 - HG6*K4)'*P*(HF6 - HG6*K4) - P <= -Q4;
            P >= eye(3);
        cvx_end       
        if strcmp(cvx_status, 'Solved');
            c = c+1;
            stab_points3(:,c) = [h,tau];
        end   
    end
end
%% Plot the stable & unstable resgions
figure(4)
hold on
for i = 1:size(points,2)
    plot(points(1,i),points(2,i),'r*');
end
for k = 1:size(stab_points3,2)
    plot(stab_points3(1,k),stab_points3(2,k),'g*');
end
xlabel('Sampling time h')
ylabel('Delay \tau')
title('Stability analysis for different combinations of h and \tau')

%% 4.4 Stability analysis from Assignemnt 1 
%Symbolic toolbox now with two variables
syms h tau

%Derive Matrix expressions and close the loop
Fx = expm(A*h);
G1 = (expm(A*(h-tau))-eye(2))/A*B;
Fu = (expm(A*(h))-eye(2))/A*B - G1;
F2 = [Fx, Fu; zeros(1,2), zeros(1)];
G2 = [G1;1];
%Extend the controller
K2 = [K, 0];
%Close the loop
A2bar_h = F2-G2*K2;

h_eval = linspace(0, 0.3, 31);
c = 0;
d = 0;
clear Abar
clear stab_points1
clear points2 
for i = 1:31
    h = h_eval(i);
    for tau = linspace(0,h,31)
        d = d+1;
        points2(:,d) = [h,tau];
        Abar = eval(A2bar_h);
        if (abs(max(eig(Abar))) <= 1);
        c = c+1;
        stab_points1(:,c) = [h,tau];
        end   
    end
end
%% 4.4 Comparison plot
figure(5)
hold on
for i = 1:size(points,2)
    plot(points(1,i),points(2,i),'r*');
end
for k = 1:size(stab_points3,2)
    plot(stab_points3(1,k),stab_points3(2,k),'g*');
end
for k = 1:size(stab_points1,2)
    plot(stab_points1(1,k),stab_points1(2,k),'y*');
end
for k = 1:size(stab_points2,2)
    plot(stab_points2(1,k),stab_points2(2,k),'b*');
end

xlabel('Sampling time h')
ylabel('Delay \tau')
title('Stability analysis for different combinations of h and \tau')

%% Question 5
%% 5.1 Quadratic event triggering condition

%Solve for Q
cvx_begin sdp 
    variable P(2,2) symmetric
    variable Q(2,2)
    subject to 
    (A-B*K)'*P+P*(A-B*K) >= -Q;
    (A-B*K)'*P+P*(A-B*K) <= -Q;
    P >= eye(2);
    Q >= eye(2);
cvx_end

%Triggering condition as a symbolic expression
syms sigma
trig_M = [(1-sigma).*Q, P*B*K; (B*K)'*P, zeros(2,2)];

%% 5.2 Simulate the closed-loop for values of sigma 
clear yout tout

sigma = 0.5;
trig_M_eval = eval(trig_M);
for k = 1:9 %Evaluate for sigma from 0.1 to 0.9
    sigma = k/10;
    trig_M_eval = eval(trig_M);
    y1 = -2.5:2.5;
    y2 = -2.5:2.5;
    [Y1,Y2] = meshgrid(y1,y2); %Different initial conditions
    for j = 1:size(Y1,1).^2;
        tstart = 0;
        tfinal = 10;
        tstop = 5;
        y0 = [Y1(j);Y2(j);0;0];
        y0 = [1;1;0;0];
        refine = 4;
        options = odeset('Events',@(t,y) events(t,y,trig_M_eval),'Refine',refine);

        tout = tstart;
        yout = y0.';
        teout = [];
        yeout = [];
        ieout = [];
        for i = 1:10000
           % Solve until the first terminal event.
           [t,y,te,ye,ie] = ode45(@(t,y) f(t,y,A,B,K),[tstart tfinal],y0,options);
           if ~ishold
              hold on
           end
           % Accumulate output.  This could be passed out as output arguments.
           nt = length(t);
           tout = [tout; t(2:nt)];
           yout = [yout; y(2:nt,:)];
           teout = [teout; te];% Events at tstart are never reported.
           yeout = [yeout; ye];
           ieout = [ieout; ie];

           if te >= tstop
               disp('Stopped')
               break
           end
           % Set the new initial conditions, with .9 attenuation.
           y0(1) = y(nt,1);
           y0(2) = y(nt,2);
           y0(3) = 0;
           y0(4) = 0;
           % A good guess of a valid first timestep is the length of the last valid
           % timestep, so use it for faster computation.  'refine' is 4 by default.
           options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
              'MaxStep',t(nt)-t(1));
           tstart = t(nt);
        end
        initial_cond(:,j) = [Y1(j);Y2(j)];
        trig(j) = size(teout,1);
    end
    avg_trig(k) = mean(trig);
end
%% Plot the results for system trajectory
figure(5)
plot(tout,yout(:,1),'r',tout,yout(:,2),'b')
hold on
plot(teout,yeout(:,1),'ro',teout,yeout(:,2),'bo')
title('Simulated closed loop response for sigma = 0.5')
legend('x1','x2');
xlabel('time(s)')
%% Plot the results for error trajectory
figure(6)
plot(tout,yout(:,3),'g',tout,yout(:,4),'m')
legend('\epsilon_1','\epsilon_2');
title('Simulated error behavior for sigma = 0.5')
xlabel('time(s)')

%% 5.3 Simulate the closed-loop for values of sigma with Lyapunov check
sigma = 0.5;
trig_M_eval = eval(trig_M);

tstart = 0;
tfinal = 300;
refine = 4;
options = odeset('Events',@(t,y) events(t,y,trig_M_eval),'Refine',refine);

y1 = -6.5:6.5;
y2 = -6.5:6.5;
[Y1,Y2] = meshgrid(y1,y2)
for j = 1:size(Y1,1).^2;
    tstart = 0;
    tfinal = 300;
    refine = 4;
    options = odeset('Events',@(t,y) events(t,y,trig_M_eval),'Refine',refine);
    y0 = [Y1(j);Y2(j);0;0];
    V_0 = y0(1:2)'*P*y0(1:2);
    tout = tstart;
    yout = y0.';
    teout = [];
    yeout = [];
    ieout = [];
    for i = 1:500
       % Solve until the first terminal event.
       [t,y,te,ye,ie] = ode23(@(t,y) f(t,y,A,B,K),[tstart tfinal],y0,options);
       if ~ishold
          hold on
       end
       % Accumulate output.  This could be passed out as output arguments.
       nt = length(t);
       tout = [tout; t(2:nt)];
       yout = [yout; y(2:nt,:)];
       teout = [teout; te];          % Events at tstart are never reported.
       yeout = [yeout; ye];
       ieout = [ieout; ie];
       if (ye(1:2)*P*ye(1:2)' <= 0.1*V_0); %Store time and triggers once condition is reached
           %disp('Decrease reached')
           count(j) = size(teout,1);
           time(j) = te;
           h_avg(j) = te/size(teout,1);
           break
       end
       % Set the new initial conditions, with .9 attenuation.
       y0(1) = y(nt,1);
       y0(2) = y(nt,2);
       y0(3) = 0;
       y0(4) = 0;
       % A good guess of a valid first timestep is the length of the last valid
       % timestep, so use it for faster computation.  'refine' is 4 by default.
       options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
          'MaxStep',t(nt)-t(1));
       tstart = t(nt);
    end
end
h_average = mean(h_avg);
%h_average = 0.08s
%% 5.4 Time-triggered conrtoller

h = 0.08;
F = expm(A*h);
G = (expm(A*h)-eye(2))*inv(A)*B;
Abar = F-G*K;
y0 = [1;1];
n = round(tout(end)/h);
y_disc = zeros(2,n);
y_disc(:,1) = y0;
t_disc = 0:0.08:tout(end);
for i = 1:n
    y_disc(:,i+1) = Abar*y_disc(:,i);
end

figure(6)
plot(tout,yout(:,1),'r',tout,yout(:,2),'b')
hold on
%plot(teout,yeout(:,1),'ro',teout,yeout(:,2),'bo')
plot(t_disc,y_disc(1,1:end-1),'ro',t_disc,y_disc(2,1:end-1),'bo');
legend('\epsilon_1','\epsilon_2');
title('Comparison of the time- and event-triggered controller for sigma = 0.5')
xlabel('time(s)')

%% --------------------------------------------------------------------------

function dydt = f(t,y,A,B,K)
dydt = [A-B*K, -B*K; -(A-B*K), B*K]*y;
end
% --------------------------------------------------------------------------

function [value,isterminal,direction] = events(t,y,trig_M_eval)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.
value = y.'*trig_M_eval*y;     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
end
