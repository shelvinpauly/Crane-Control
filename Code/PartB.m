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
A = simplify((subs(Ja , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])))
B = simplify((subs(Jb , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])))
C = simplify((subs(Jc , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])))
D = simplify((subs(Jd , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])))