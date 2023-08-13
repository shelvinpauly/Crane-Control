
%% Define the state space equation
syms F m1 m2 t1 td1 L2 L1 t2 td2 g M xd x
f1 = xd;
f3 = td1;
f5 = td2;
f2 = (F - m1*L1*sin(t1)*td1^2 - m2*L2*sin(t2)*td2^2 - m1*g*sin(t1)*cos(t1) - m2*g*sin(t2)*cos(t2))/(M + m1 + m2 - m1*cos(t1)^2 - m2*cos(t2)^2);
f4 = cos(t1)*f2/L1 - g*sin(t1)/L1;
f6 = cos(t2)*f2/L2 - g*sin(t2)/L2;
tspan = 0:0.1:100;
q0 = [5 0 deg2rad(0) 0 deg2rad(0) 0];


%% Find the Jacobian to Linearize the system.
%% A, B, C and D matrices are found using jacobian
Ja = jacobian ([f1, f2, f3, f4, f5, f6],[x, xd, t1, td1, t2, td2]);
Jb = jacobian ([f1, f2, f3, f4, f5, f6],[F]);
Jc = jacobian ([x, t1, t2],[x, xd, t1, td1, t2, td2]);
Jd = jacobian ([x, t1, t2], [F]);



%% Find the A, B, C and D matrix around equilibrium point which is (0,0)
A = simplify((subs(Ja , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])))
B = simplify((subs(Jb , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
C = simplify((subs(Jc , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
D = simplify((subs(Jd , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));

%% Linearized Model
A0 = double(subs(A,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
B0 = double(subs(B,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
C0 = double(subs(C,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
D0 = double(subs(D,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));


%%check the observality of the system when output is x
Jc2 = jacobian ([x, 0, 0],[x, xd, t1, td1, t2, td2]);
C2 = double(simplify((subs(Jc2 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))))

sys2 = ss(A0, B0, C2, D0);

Q = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
R = 0.3;
[K1,S,P] = lqr(A0, B0, Q, R);

sys = ss(A0-B0*K1,B0,C2,D0);

%% Kalman Estimator Design
Bd = 0.5*eye(6);                %Process Noise
Vn = 0.05;                      %Measurement Noise
[L2,P,E] = lqe(A0,Bd,C2,Bd,Vn*eye(3)); %Considering vector output: x(t)
Ac1 = A0-(L2*C2);
e_sys1 = ss(Ac1,[B0 L2],C2,0);

%% Non-linear Model LQG Response
[t,q1] = ode45(@(t,q)nonLinear(t,q,-K1*q,L2),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variable')
xlabel('time (sec)')
title('Non-Linear System LQG for output vector: x(t)')
legend('x')
hold off

function dQ = nonLinear(t,y,F,Lue1)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    dx = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y1 = [x; 0; 0];
    c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
    sum = Lue1*(y1-c1*y);
    dQ(1) = dx + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end