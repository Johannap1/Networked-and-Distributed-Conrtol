%% Initially given
%Student number: 4677498
clear all
clc
a = 4;
b = 7;
c = 8;
A = [a b+0.5; 0 -c];
B = [0;1];

%% -------------- Question 1 ---------------------
% Constant sampling interval h, no dealy tau 

%% 1. Place poles at -2, -3 for continuous system
p = [-2 -3];
K = place(A,B,p);

%% 2. Construct the discrete-time model

%Symbolic toolbox to get an expression as a function of h
syms h 

%Compute matrices F & G
F = expm(A*h);
G = (expm(A*h)-eye(2))*inv(A)*B;

%close the loop
Abar_h = F-G*K;

%% 3. Stability as a function of sampling interval h 

%Initialize arrays
EV = zeros(2,100);
max_eig = zeros(1,100);
c = 0;

%Compute maximum Eigenvalue for different values of h
for h = linspace(0,1,100);
    c = c+1;
    Abar = eval(Abar_h);
    EV(:,c) = eig(Abar);
    max_eig(c) = abs(max(eig(Abar)));
end

%Plot the maximum Eigenvalue as a function of h
figure(1)
h = linspace(0,1,100);
plot(h, max_eig);
hold on
plot(h,ones(1,100),'r--');
xlabel('Sampling Interval h')
ylabel('|\lambda_{max}|')
title('Norm of the maximum Eigenvalue')

%Plot the Evolution of the Eigenvalues
figure(2)
plot(EV(1,:),'r');
hold on
plot(EV(2,:),'g');
xlabel('Real(\lambda)')
ylabel('Imag(\lambda)')
title('Eigenvalues')

%% --------------- Question 2 -------------------
% Different comninations of sampling intervals and time delays
%% 1. Construct the model

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

%% 2. Study Stability

%Compute stable and unstable points for different combinations of h and tau
h_eval = linspace(0, 0.3, 31);
c = 0;
d = 0;
clear Abar
clear stab_points2
clear points2 
for i = 1:31
    h = h_eval(i);
    for tau = linspace(0,h,31)
        d = d+1;
        points2(:,d) = [h,tau];
        Abar = eval(A2bar_h);
        if (abs(max(eig(Abar))) <= 1);
        c = c+1;
        stab_points2(:,c) = [h,tau];
        end   
    end
end

% Plot the stable vs unstable points
figure(3)
hold on
for i = 1:size(points2,2)
    plot(points2(1,i),points2(2,i),'r*');
end
for k = 1:size(stab_points2,2)
    plot(stab_points2(1,k),stab_points2(2,k),'g*');
end
set(gca, 'Fontsize', 15)
xlabel('Sampling time h','Fontsize', 20)
ylabel('Delay \tau','Fontsize', 20)
title('Stability analysis for different combinations of h and \tau','Fontsize', 20)

%% 3. Controller for robustness against delays

%chose specific point at which the closed loop is evaluated
h = 0.2;
tau = 0.2;
A23 = eval(F2);
B23 = eval(G2);
%place the poles
p = [0.91,0.92,0.9];
K23 = place(A23,B23,p);
%close the loop with the new controller
A23bar_h = F2-G2*K23;

%Study stability again 
h_eval = linspace(0, 0.5, 51);
c = 0;
d = 0;
clear Abar
clear stab_points23
clear points23 
for i = 1:51
    h = h_eval(i);
    for tau = linspace(0,h,51)
        d = d+1;
        points23(:,d) = [h,tau];
        Abar = eval(A23bar_h);
        if (abs(max(eig(Abar))) <= 1);
        c = c+1;
        stab_points23(:,c) = [h,tau];
        end   
    end
end

%plot the stable and unstable points
figure(5)
hold on
for i = 1:size(points23,2)
    plot(points23(1,i),points23(2,i),'r*');
end
for k = 1:size(stab_points23,2)
    plot(stab_points23(1,k),stab_points23(2,k),'g*');
end
set(gca, 'Fontsize', 15)
xlabel('Sampling time h','Fontsize', 20)
ylabel('Delay \tau','Fontsize', 20)
title('Stability analysis for different combinations of h and \tau','Fontsize', 20)

%% ------------ Question 3 -------------------- %%
%% 1. Construct the model
%Symbolic expression
syms h tau

%Construct the system matrices
Fx = expm(A*h);
Fu1 = (expm(A*(2*h-tau))-eye(2))/A*B;
Fu2 = (expm(A*(h))-eye(2))/A*B - Fu1;
F3 = [Fx, Fu1, Fu2; zeros(1,4); zeros(1,2), 1, 0];
G3 = [zeros(2,1);1;0];
%Extend the controller
K3 = [K, 0, 0];
%Close the loop
A3bar_h = F3-G3*K3;

%% 2. Study Stability
h_eval = linspace(0, 0.3, 31);
c = 0;
d = 0;
clear Abar
clear stab_points3
clear points3
for i = 1:31
    h = h_eval(i);
    for tau = linspace(h,2*h,31)
        d = d+1;
        points3(:,d) = [h,tau];
        Abar = eval(A3bar_h);
        if (abs(max(eig(Abar))) <= 1);
        c = c+1;
        stab_points3(:,c) = [h,tau];
        end   
    end
end
% Plot the stable and unstable points on top of the results from question 2
figure(3)
hold on
for i = 1:size(points2,2)
    plot(points2(1,i),points2(2,i),'r*');
end
for k = 1:size(stab_points2,2)
    plot(stab_points2(1,k),stab_points2(2,k),'g*');
end
for i = 1:size(points3,2)
    plot(points3(1,i),points3(2,i),'r*');
end
for k = 1:size(stab_points3,2)
    plot(stab_points3(1,k),stab_points3(2,k),'g*');
end
set(gca, 'Fontsize', 15)
xlabel('Sampling time h','Fontsize', 20)
ylabel('Delay \tau','Fontsize', 20)
title('Stability analysis for different combinations of h and \tau','Fontsize', 20)

%% 3. Controller for robustness against delays
%Approach 1: Pole placcement
tau = 0.1;
h = 0.2;
A33 = eval(F3);
B33 = G3;
%p = [0.95,0.9,0.75,0.7];
%K33 = place(A33,B33,p);

%Approach 2: Dynamic conrtoller form 2.3
K33 = [K23 0];
A33bar_h = F3-G3*K33;

%% Test again for stability
h_eval = linspace(0, 0.5, 51);
c = 0;
d = 0;
clear Abar
clear stab_points33
clear points33
for i = 1:51
    h = h_eval(i);
    for tau = linspace(h,2*h,51)
        d = d+1;
        points33(:,d) = [h,tau];
        Abar = eval(A33bar_h);
        if (abs(max(eig(Abar))) <= 1);
        c = c+1;
        stab_points33(:,c) = [h,tau];
        end   
    end
end

%plot the results
figure(8)
hold on
for i = 1:size(points23,2)
    plot(points23(1,i),points23(2,i),'r*');
end
for k = 1:size(stab_points23,2)
    plot(stab_points23(1,k),stab_points23(2,k),'g*');
end
for i = 1:size(points33,2)
    plot(points33(1,i),points33(2,i),'r*');
end
for k = 1:size(stab_points33,2)
    plot(stab_points33(1,k),stab_points33(2,k),'g*');
end
set(gca, 'Fontsize', 15)
xlabel('Sampling time h','Fontsize', 20)
ylabel('Delay \tau','Fontsize', 20)
title('Stability analysis for different combinations of h and \tau','Fontsize', 20)

%% ----------------- Question 4 ---------------- 
%% 1. Test for stability
%Evaluate for increasing values of h
h_eval = linspace(0, 0.3, 31);
d = 0;
clear stab_points4
clear points4
Q = eye(4);
for i = 1:31
    h = h_eval(i);
    tau = 0.2*h;
    A1 = [eval(A2bar_h),zeros(3,1);0 0 1 0];
    tau = 0.5*h;
    A2 = [eval(A2bar_h),zeros(3,1); 0 0 1 0];
    tau = h;
    A3 = eval(A3bar_h);
    tau = 1.5*h;
    A4 = eval(A3bar_h);
    cvx_begin sdp 
    variable P(4,4) symmetric 
    subject to 
        A1'*P*A1 - P <= -Q;
        A2'*P*A2 - P <= -Q;
        A3'*P*A3 - P <= -Q;
        A4'*P*A4 - P <= -Q;
        P >= eye(4);
    cvx_end
    if (strcmp(cvx_status, 'Solved'))
        d = d+1;
        stab_points4(d) = h; %Store h in vector if system is stable
    end
end

%% 2. Design a controller

h_eval = linspace(0.01, 0.3, 30);
Q = eye(8);
for i = 1:31
    h = h_eval(i);
    tau = 0.2*h;
    F41 = [eval(F2),zeros(3,1);0 0 1 0];
    G41 = [eval(G2);0];
    tau = 0.5*h;
    F42 = [eval(F2),zeros(3,1);0 0 1 0];
    G42 = [eval(G2);0];
    tau = h;
    F43 = eval(F3);
    G43 = G3;
    tau = 1.5*h;
    F44 = eval(F3);
    G44 = G3;
    cvx_begin sdp 
    variable M(4,4)
    variable Y(1,4)
    subject to 
        [-M M'*F41'-Y'*G41';F41*M-G41*Y -M] <= -Q;
        [-M M'*F42'-Y'*G42';F42*M-G42*Y -M] <= -Q;
        [-M M'*F43'-Y'*G43';F43*M-G43*Y -M] <= -Q;
        [-M M'*F44'-Y'*G44';F44*M-G44*Y -M] <= -Q;
    cvx_end
    P = inv(M);
    if (~strcmp(cvx_status, 'Solved'))
        h
        break
    end
end

%% Choose the last value of h for which the system was stable and use that controller
h = 0.16;
tau = 0.2*h;
F41 = [eval(F2),zeros(3,1);0 0 1 0];
G41 = [eval(G2);0];
tau = 0.5*h;
F42 = [eval(F2),zeros(3,1);0 0 1 0];
G42 = [eval(G2);0];
tau = h;
F43 = eval(F3);
G43 = G3;
tau = 1.5*h;
F44 = eval(F3);
G44 = G3;
Q = eye(8);

cvx_begin sdp 
    variable M(4,4) semidefinite
    variable Y(1,4)
    subject to 
        [-M M'*F41'-Y'*G41';F41*M-G41*Y -M] <= -Q;
        [-M M'*F42'-Y'*G42';F42*M-G42*Y -M] <= -Q;
        [-M M'*F43'-Y'*G43';F43*M-G43*Y -M] <= -Q;
        [-M M'*F44'-Y'*G44';F44*M-G44*Y -M] <= -Q;
cvx_end
P = inv(M);
K4 = Y/M;
%% Check for stability again

h_eval = linspace(0, 0.3, 31);
d = 0;
clear stab_points4
clear points4
Q = eye(4);
for i = 1:31
    h = h_eval(i);
    tau = 0.2*h;
    F41 = [eval(F2),zeros(3,1);0 0 1 0];
    G41 = [eval(G2);0];
    A1 = F41 - G41*K4;
    tau = 0.5*h;
    F42 = [eval(F2),zeros(3,1);0 0 1 0];
    G42 = [eval(G2);0];
    A2 = F42 - G42*K4;
    tau = h;
    F43 = eval(F3);
    G43 = G3;
    A3 = F43 - G43*K4;
    tau = 1.5*h;
    F44 = eval(F3);
    G44 = G3;
    A4 = F44 - G44*K4;
    cvx_begin sdp 
    variable P(4,4) symmetric 
    %variable Q(4,4) symmetric 
    subject to 
        A1'*P*A1 - P <= -Q;
        A2'*P*A2 - P <= -Q;
        A3'*P*A3 - P <= -Q;
        A4'*P*A4 - P <= -Q;
        P >= eye(4);
        %Q >= eye(4);
    cvx_end
    if (strcmp(cvx_status, 'Solved'))
        d = d+1;
        stab_points4(d) = h;
    end
end
%% 3. Check stability for periodic system

c = 0;
clear Abar

for h = linspace(0, 0.3, 31);
    tau = 0.2*h;
    A1 = [eval(A2bar_h),zeros(3,1);0 0 1 0];
    tau = 0.5*h;
    A2 = [eval(A2bar_h),zeros(3,1); 0 0 1 0];
    tau = h;
    A3 = eval(A3bar_h);
    c = c+1;
    Abar = A2*A3*A1;
    max_eig(c) = abs(max(eig(Abar)));  
end

figure(7)
h = linspace(0, 0.3, 31);
plot(h, max_eig);
hold on
plot(h,ones(1,31),'r--');
xlabel('Sampling Interval h')
ylabel('|\lambda_{max}|')
title('Norm of the maximum Eigenvalue')

%% 4. Design a controller for the periodic system

% Same approach as before now with only three system matrices
h_eval = linspace(0.01, 0.3, 30);
Q = eye(8);
for i = 1:30
    h = h_eval(i);
    tau = 0.2*h;
    F41 = [eval(F2),zeros(3,1);0 0 1 0];
    G41 = [eval(G2);0];
    tau = 0.5*h;
    F42 = [eval(F2),zeros(3,1);0 0 1 0];
    G42 = [eval(G2);0];
    tau = h;
    F43 = eval(F3);
    G43 = G3;
    tau = 1.5*h;
    F44 = eval(F3);
    G44 = G3;
    cvx_begin sdp 
    variable M(4,4)
    variable Y(1,4)
    subject to 
        [-M M'*F41'-Y'*G41';F41*M-G41*Y -M] <= -Q;
        [-M M'*F42'-Y'*G42';F42*M-G42*Y -M] <= -Q;
        [-M M'*F43'-Y'*G43';F43*M-G43*Y -M] <= -Q;
        %[-M M'*F44'-Y'*G44';F44*M-G44*Y -M] <= -Q;
    cvx_end
    P = inv(M);
    if (~strcmp(cvx_status, 'Solved'))
        h
        break
    end
end

%% Find the controller for the found h
h = 0.24;
tau = 0.2*h;
F41 = [eval(F2),zeros(3,1);0 0 1 0];
G41 = [eval(G2);0];
tau = 0.5*h;
F42 = [eval(F2),zeros(3,1);0 0 1 0];
G42 = [eval(G2);0];
tau = h;
F43 = eval(F3);
G43 = G3;

cvx_begin sdp 
    variable M(4,4)
    variable Y(1,4)
    subject to 
        [-M M'*F41'-Y'*G41';F41*M-G41*Y -M] <= -Q;
        [-M M'*F42'-Y'*G42';F42*M-G42*Y -M] <= -Q;
        [-M M'*F43'-Y'*G43';F43*M-G43*Y -M] <= -Q;
        cvx_end
P = inv(M);
K42 = Y/M;

%% Reevaluate stability for both the new controller K42 as well as the one from Question 2 (K4)

% Old controller
c = 0;
clear Abar
clear max_eig
for h = linspace(0, 1, 101);
    tau = 0.2*h;
    F41 = [eval(F2),zeros(3,1);0 0 1 0];
    G41 = [eval(G2);0];
    A1 = F41 - G41*K4;
    tau = 0.5*h;
    F42 = [eval(F2),zeros(3,1);0 0 1 0];
    G42 = [eval(G2);0];
    A2 = F42 - G42*K4;
    tau = h;
    F43 = eval(F3);
    G43 = G3;
    A3 = F43 - G43*K4;
    c = c+1;
    Abar =  A2*A3*A1;
    max_eig(c) = abs(max(eig(Abar)));  
end

% New controller
c = 0;
clear Abar
clear max_eig_2
for h = linspace(0, 1, 101);
    tau = 0.2*h;
    F41 = [eval(F2),zeros(3,1);0 0 1 0];
    G41 = [eval(G2);0];
    A1 = F41 - G41*K42;
    tau = 0.5*h;
    F42 = [eval(F2),zeros(3,1);0 0 1 0];
    G42 = [eval(G2);0];
    A2 = F42 - G42*K42;
    tau = h;
    F43 = eval(F3);
    G43 = G3;
    A3 = F43 - G43*K42;
    c = c+1;
    Abar =  A2*A3*A1;
    max_eig_2(c) = abs(max(eig(Abar)));  
end

% plot the results
figure(8)
h = linspace(0, 1, 101);
plot(h, max_eig);
hold on
plot(h, max_eig_2);
ylim([-0.5,1.5])
plot(h,ones(1,101),'r--');
xlabel('Sampling Interval h')
ylabel('|\lambda_{max}|')
title('Norm of the maximum Eigenvalue')
legend('K_{4.2}','K_{4.4}')
