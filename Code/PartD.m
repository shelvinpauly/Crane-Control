%% Define the state space equation
syms F m1 m2 t1 td1 L2 L1 t2 td2 g M xd x
f1 = xd;
f3 = td1;
f5 = td2;
f2 = (F - m1*L1*sin(t1)*td1^2 - m2*L2*sin(t2)*td2^2 - m1*g*sin(t1)*cos(t1) - m2*g*sin(t2)*cos(t2))/(M + m1 + m2 - m1*cos(t1)^2 - m2*cos(t2)^2);
f4 = cos(t1)*f2/L1 - g*sin(t1)/L1;
f6 = cos(t2)*f2/L2 - g*sin(t2)/L2;


%% Find the Jacobian to Linearize the system.
%% A, B, C and D matrices are found using jacobian
Ja = jacobian ([f1, f2, f3, f4, f5, f6],[x, xd, t1, td1, t2, td2]);
Jb = jacobian ([f1, f2, f3, f4, f5, f6],[F]);
Jc = jacobian ([x, t1, t2],[x, xd, t1, td1, t2, td2]);
Jd = jacobian ([x, t1, t2], [F]);



%% Find the A, B, C and D matrix around equilibrium point which is (0,0)
A = simplify((subs(Ja , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
B = simplify((subs(Jb , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
C = simplify((subs(Jc , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
D = simplify((subs(Jd , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));



%% Find the rank of the controllability matrix of the linearized system around (0,0)
%% A, B, C and D matrix are symbolically represented during this rank calculation
E = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];
rank(E);


%% check that the system is controllable when M = 1000, m1 = m2 = 100, L1 = 20, L2 = 10
E1 = subs(E, [M m1 m2 L1 L2], [1000 100 100 20 10]);
disp("Rank of the controllability matrix: ")
rank(E1)

%% LQR design

%%LQR without gain tuning
Q = eye(6,6);
R = 1;
A0 = double(subs(A,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
B0 = double(subs(B,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
C0 = double(subs(C,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
D0 = double(subs(D,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
K1 = lqr(A0, B0, Q, R);

sys = ss((A0-B0*K1),B0, C0, D0);
figure(1)
step(sys)

%% LQR after gain tuning
Q(4,4) = 1000;
Q(6,6) = 1000;
R = 0.0001;
A1 = double(subs(A,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
B1 = double(subs(B,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
C1 = double(subs(C,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
D1 = double(subs(D,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
K2 = lqr(A1, B1, Q, R);

tspan = 0:0.1:100;
X0 = [5; 0; deg2rad(10); 0; deg2rad(10); 0];  %% intial condition of the system
Xw = [0; 0; 0; 0; 0; 0];  %% Final condition of the system
u = @(X) -K2*(X - Xw);     %% Control input
[t,X] = ode45(@(t,X)system(X, A1, B1, u(X)), tspan, X0);

title('LQR on the linear controller')
figure(2)
subplot(3,1,1);
plot(t, X(:,1))
xlabel('time(s)');
ylabel('x');
subplot(3,1,2);
plot(t, X(:,2))
xlabel('time(s)');
ylabel('theta1');
subplot(3,1,3);
plot(t, X(:,3))
xlabel('time(s)')
ylabel('theta2')

% %% Non_Linear Model Controlled Response
% tspan = 0:0.1:100;
% X0 = [5; 0; deg2rad(10); 0; deg2rad(10); 0];  %% intial condition of the system
% Xw = [0; 0; 0; 0; 0; 0];  %% Final condition of the system
% [t,q2] = ode45(@(t,q)LQRonNonLinear(t,q,-K2*q),tspan,X0);
% figure(2);
% hold on
% plot(t,X1(:,1))
% plot(t,X1(:,3))
% plot(t,X1(:,5))
% ylabel('state variables')
% xlabel('time (sec)')
% title('Non-Linear system using LQR controller')
% legend('x','theta1','theta2')

%% Lyapunov indirect method to certify stability of the close loop system
%% Eigen value of the close loop system is calculated
%% since all of the eigen value has a negative real part, system is stable 
Ac = A1 - B1*K2;
disp("eigen value of the closed loop matrix:")
ei = eig(Ac)
disp("since all of the eigen value has negative real part, system is controllable")

function dX = system(X, A1, B1, u)
    dX = A1*X + B1*u;
end

% function dQ = LQRonNonLinear(t,y,F)
%     m1 = 100; m2 = 100; M=1000; L1 = 20; L2 = 10; g = 9.81;
%     x = y(1);
%     dx = y(2);
%     t1 = y(3);
%     td1 = y(4);
%     t2 = y(5);
%     td2 = y(6);
%     dQ = zeros(6:1);
%     dQ(1) = dx;
%     dQ(2) = (F - m1*L1*sin(t1)*td1^2 - m2*L2*sin(t2)*td2^2 - m1*g*sin(t1)*cos(t1) - m2*g*sin(t2)*cos(t2))/(M + m1 + m2 - m1*cos(t1)^2 - m2*cos(t2)^2);
%     dQ(3) = td1;
%     dQ(4) = cos(t1)*dQ(2)/L1 - g*sin(t1)/L1;
%     dQ(5) = td2;
%     dQ(6) = cos(t2)*dQ(2)/L2 - g*sin(t2)/L2;
% end
    