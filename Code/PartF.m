%% Define the state space equation
syms F m1 m2 t1 td1 L2 L1 t2 td2 g M xd x
f1 = xd;
f3 = td1;
f5 = td2;
f2 = (F - m1*L1*sin(t1)*td1^2 - m2*L2*sin(t2)*td2^2 - m1*g*sin(t1)*cos(t1) - m2*g*sin(t2)*cos(t2))/(M + m1 + m2 - m1*cos(t1)^2 - m2*cos(t2)^2);
f4 = cos(t1)*f2/L1 - g*sin(t1)/L1;
f6 = cos(t2)*f2/L2 - g*sin(t2)/L2;
X0 = [5; 0; deg2rad(5); 0; deg2rad(10); 0]

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


%% 
Q = eye(6,6);
R = 1;
A0 = double(subs(A,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
B0 = double(subs(B,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
C0 = double(subs(C,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
D0 = double(subs(D,[M m1 m2 L1 L2 g], [1000 100 100 20 10 9.8]));
K1 = lqr(A0, B0, Q, R);

%%check the observality of the system when output is x
Jc2 = jacobian ([x, 0, 0],[x, xd, t1, td1, t2, td2]);
C2 = double(simplify((subs(Jc2 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))))

%%check the observality of the system when output is x and t2
Jc3 = jacobian ([x, 0, t2],[x, xd, t1, td1, t2, td2]);
C3 = double(simplify((subs(Jc3 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))));

%%check the observality of the system when output is x t1 and t2
Jc4 = jacobian ([x, t1, t2],[x, xd, t1, td1, t2, td2]);
C4 = double(simplify((subs(Jc4 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))));

sys2 = ss(A0, B0, C2, D0);
sys3 = ss(A0, B0, C3, D0);
sys4 = ss(A0, B0, C4, D0);


%% kalman estimation
Bd = 0.1*eye(6,6);
Vt = 0.01*eye(3,3);
[L2,P, E] = lqe(A0, Bd, C2, Bd,Vt);
[L3,P, E] = lqe(A0, Bd, C3, Bd,Vt);
[L4,P, E] = lqe(A0, Bd, C4, Bd,Vt);

Ac2 = A0 - L2*C2;
Ac3 = A0 - L3*C3;
Ac4 = A0 - L4*C4;

e_sys2 = ss(Ac2, [B0 L2], C2, 0);
e_sys3 = ss(Ac3, [B0 L3], C3, 0);
e_sys4 = ss(Ac4, [B0 L4], C4, 0);
%% step input response
tspan = 0:0.1:100;
unitStep = 0*tspan;
unitStep(200:length(tspan)) = 1;

[Y2,t] = lsim(sys2,unitStep, tspan);
[X2,t] = lsim(e_sys2,[unitStep;Y2'],tspan);

[Y3,t] = lsim(sys3,unitStep, tspan);
[X3,t] = lsim(e_sys3,[unitStep;Y3'],tspan);

[Y4,t] = lsim(sys4,unitStep, tspan);
[X4,t] = lsim(e_sys2,[unitStep;Y4'],tspan);

figure();
hold on
plot(t,Y2(:,1),'r','Linewidth',2)
plot(t,X2(:,1),'k--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','Estimated x(t)')
title('Response for output vector at step input: (x(t)')
hold off

figure();
hold on
plot(t,Y3(:,1),'r','Linewidth',2)
plot(t,Y3(:,3),'b','Linewidth',2)
plot(t,X3(:,1),'k--','Linewidth',1)
plot(t,X3(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_2(t)','Estimated x(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_2(t))')
hold off

figure();
hold on
plot(t,Y4(:,1),'r','Linewidth',2)
plot(t,Y4(:,2),'g','Linewidth',2)
plot(t,Y4(:,3),'b','Linewidth',2)
plot(t,Y4(:,1),'k--','Linewidth',1)
plot(t,X4(:,2),'r--','Linewidth',1)
plot(t,X4(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('time(sec)')
legend('x(t)','theta_1(t)','theta_2(t)','Estimated x(t)','Estimated theta_1(t)','Estimated theta_2(t)')
title('Response for output vector at step input: (x(t),theta_1(t),theta_2(t))')
hold off

[t,x2] = ode45(@(t,x)linear2(t, x  ,L2, A0, B0, C2),tspan,X0);
figure();
hold on
plot(t,x2(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Linear system Observer when output vector: x(t)')
legend('x')
hold off

[t,x3] = ode45(@(t,x)linear3(t, x, L3, A0, B0, C3),tspan,X0);
figure();
hold on
plot(t,x3(:,1))
plot(t,x3(:,5))
ylabel('state variables')
xlabel('time (s)')
title('Linear system Observer when output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off

[t,x4] = ode45(@(t,x)linear4(t, x, L4, A0, B0, C4),tspan,X0);
figure();
hold on
plot(t,x4(:,1))
plot(t,x4(:,3))
plot(t,x4(:,5))
ylabel('state variables')
xlabel('time (s)')
title('Linear system Observer when output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

%% Non-linear Model Observer Response
[t,q1] = ode45(@(t,q)nonLinear2(t,q,1,L2),tspan,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)nonLinear3(t,q,1,L3),tspan,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_2(t))')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)nonLinear4(t,q,1,L4),tspan,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear System Observer for output vector: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

function dX = linear2(t,x,L,A,B,C)
    y = [x(1); 0; 0];
    K = 1;
    dX = (A+B*K)*x + L*(y - C*x);
end

function dX = linear3(t,x,L,A,B,C)
    y = [x(1); 0; x(5)];
    K = 1;
    dX = (A+B*K)*x + L*(y - C*x);
end
function dX = linear4(t,x,L,A,B,C)
    y = [x(1); x(3); x(5)];
    K = 1;
    dX = (A+B*K)*x + L*(y - C*x);
end

function dQ = nonLinear2(t,y,F,L2)
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
    sum = L2*(y1-c1*y);
    dQ(1) = dx + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

function dQ = nonLinear3(t,y,F,L3)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    dx = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y3 = [x; 0; t2];
    c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
    sum = L3*(y3-c3*y);
    dQ(1) = dx + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

function dQ = nonLinear4(t,y,F,Lue4)
    m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
    x = y(1);
    dx = y(2);
    t1 = y(3);
    dt1 = y(4);
    t2 = y(5);
    dt2 = y(6);
    dQ=zeros(6,1);
    y4 = [x; t1; t2];
    c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
    sum = Lue4*(y4-c4*y);
    dQ(1) = dx + sum(1);
    dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
    dQ(3) = dt1+sum(3);
    dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
    dQ(5) = dt2 + sum(5);
    dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end

