V_in=9;
L=.3e-3;
C=.07e-3;
R=10;

D_ss = 0.6;
VO_ss = V_in/(1-D_ss);
iL_ss = VO_ss/(R*(1-D_ss));



A = [0 -(1-D_ss)/L;
    (1-D_ss)/C -1/(R*C)];
B = [VO_ss/L;
    -iL_ss/C];