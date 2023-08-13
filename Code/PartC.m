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
disp("Rank of the controllability matrix of the linearized system: ")
rank(E)



%%Check the condition when the system is not controllable
%% when and L1 = L2, controllability matrix loose its rank, hence system is not controllable
%% Other than that, system is controllable for all other the conditions
E1 = subs(E, [m1 m2 L1 L2], [100 160 50 50]);
disp("Rank of the controllability matrix of the linearized system when length are equal: ")
rank(E1)
disp("Since rank is less than 6, system is not controllable. ")