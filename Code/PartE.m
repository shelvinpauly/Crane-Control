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

%%check the observality of the system when output is x
Jc2 = jacobian ([x, 0, 0],[x, xd, t1, td1, t2, td2]);
C2 = double(simplify((subs(Jc2 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))))
O1 = [C2 ; C2*A ; C2*A^2 ; C2*A^3 ; C2*A^4 ; C2*A^5];
disp("Rank of the observability matrix when output vector is x : ")
rank(O1)
disp("since the rank is 6, system is observable.")

%%check the observality of the system when output is t1 and t2
Jc3 = jacobian ([0, t1, t2],[x, xd, t1, td1, t2, td2]);
C3 = simplify((subs(Jc3 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0])));
O2 = [C3 ; C3*A ; C3*A^2 ; C3*A^3 ; C3*A^4 ; C3*A^5];
disp("Rank of the observability matrix when output vector is theta1 and theta2")
rank(O2)
disp("since the rank is less than 6, system is not observable.")

%%check the observality of the system when output is x and t2
Jc3 = jacobian ([x, 0, t2],[x, xd, t1, td1, t2, td2]);
C3 = double(simplify((subs(Jc3 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))));
O3 = [C3 ; C3*A ; C3*A^2 ; C3*A^3 ; C3*A^4 ; C3*A^5];
disp("Rank of the observability matrix when output vector is x and theta2")
rank(O3)
disp("since the rank is 6, system is observable.")

%%check the observality of the system when output is x t1 and t2
Jc4 = jacobian ([x, t1, t2],[x, xd, t1, td1, t2, td2]);
C4 = double(simplify((subs(Jc4 , [t1 td1 t2 td2 xd x], [0 0 0 0 0 0]))));
O4 = [C4 ; C4*A ; C4*A^2 ; C4*A^3 ; C4*A^4 ; C4*A^5];
disp("Rank of the observability matrix when output vector is x theta1 and theta2")
rank(O4)
disp("since the rank is 6, system is observable.")