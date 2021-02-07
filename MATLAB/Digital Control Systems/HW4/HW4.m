
%Problem 1
G_s = tf([1],[1 1 0]);
G_z = c2d(G_s,0.02)
p_Gs = pole(G_s)
z_Gs = zero(G_s)
p_Gz = pole(G_z)
z_Gz = zero(G_z)

%Problem 2
A = [(-5/2) 1;-175 0]
B = [0;80]
C = [1 0]
sI = [s 0; 0 s]
Gs = C*((sI-A)^-1)*B
G_s = tf([160],[2 5 350])
G_z = c2d(G_s,0.00252)

%Problem 3
G_s = tf([3600],[1 84.853 3600])
G_z = c2d(G_s, 0.0035)
p_Gs = pole(G_s)
z_Gs = zero(G_s)
p_Gz = pole(G_z)
z_Gz = zero(G_z)

