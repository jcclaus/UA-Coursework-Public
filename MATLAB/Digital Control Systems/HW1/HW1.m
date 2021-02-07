% John Claus
% 01/15/2017
% ECE 542
% HW#1 MATLAB Assignment


%Part a - First Order System Response
ta=0:0.01:10.0;
nca=1;
dca=[0,0,1];
gca=tf(nca,dca);
figure(1)
step(gca,ta)
title('Part (a): First Order System Response')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part b - Critically Damped System Response
tb=0:0.01:10.0;
ncb=1;
dcb=[1,2,1];
gcb=tf(ncb,dcb);
figure(2)
step(gcb,tb)
title('Part (b): Critically Damped System Response')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part c - Overdamped System Response
tc=0:0.01:10.0;
ncc=4;
dcc=[1,5,4];
gcc=tf(ncc,dcc);
figure(3)
step(gcc,tc)
title('Part (c): Overdamped System Response')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Part d - No Finite Zeroes w/ Varying Zeta
td=0:0.01:50.0;
ncd=1;

%Subset i - Zeta=0.707
zeta1=0.707;
dcd1=[1,(2*zeta1),1];
gcd1=tf(ncd,dcd1);
figure(4)
step(gcd1,td)
title('Part (d)(i): No Finite Zeroes (zeta=0.707)')

%Subset ii - Zeta=0.45
zeta2=0.45;
dcd2=[1,(2*zeta2),1];
gcd2=tf(ncd,dcd2);
figure(5)
step(gcd2,td)
title('Part (d)(ii): No Finite Zeroes (zeta=0.45)')

%Subset iii - Zeta=0.1
zeta3=0.1;
dcd3=[1,(2*zeta3),1];
gcd3=tf(ncd,dcd3);
figure(6)
step(gcd3,td)
title('Part (d)(iii): No Finite Zeroes (zeta=0.1)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part e - Different Zero Locations 
te=0:0.01:10.0;
dce=[1,(2*0.707), 1];

%Subset i - z=5
z1=5;
nce1=[(1/z1),1];
gce1=tf(nce1,dce);
figure(7)
step(gce1,te)
title('Part (e)(i): Different Zero Locations (z=5)')

%Subset ii - z=1
z2=1;
nce2=[(1/z2),1];
gce2=tf(nce2,dce);
figure(8)
step(gce2,te)
title('Part (e)(ii): Different Zero Locations (z=1)')

%Subset iii - z=0.4
z3=0.4;
nce3=[(1/z3),1];
gce3=tf(nce3,dce);
figure(9)
step(gce3,te)
title('Part (e)(iii): Different Zero Locations (z=0.4)')

%Subset iv - z=-1
z4=-1;
nce4=[(1/z4),1];
gce4=tf(nce4,dce);
figure(10)
step(gce4,te)
title('Part (e)(iv): Different Zero Locations (z=-1)')
