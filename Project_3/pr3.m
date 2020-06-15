G_CL  = tf([5 0 0],[-25 40 2000 0 0]);
C2tf = tf(C);
Cr = C2tf(1);
Cy = C2tf(2);
T  = Cr*feedback(G_CL,Cy,1);
step(T);
grid on;
C