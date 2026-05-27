%==========================================================================
% MAUB Computation via Theorem 1
% System: x(k+1) = B*x(k) + Bd*x(k-d(k)),  d_l <= d(k) <= d_u
% Toolbox: YALMIP + sdpt3
%==========================================================================
close;
clear;
clc;


%--------------------------------------------------------------------------
% System Selection
%--------------------------------------------------------------------------
example = 1;   % Set to 1 or 2
switch example
    case 1
        B  = [1.00,  0.01; -0.10, 0.99];
        Bd = [0.003, 0.001; 0.010, 0.005];
    case 2
        B  = [0.80, 0.00;  0.05, 0.90];
        Bd = [-0.1, 0.0; -0.2, -0.1];
end
n=2;

% Initialize basis matrices e0-e18
e0=zeros(n,18*n);
basisSize = 18; % 18n-dimensional basis
e = cell(1, basisSize);
for i = 1:basisSize
    e{i} = [zeros(n, (i-1)*n) eye(n) zeros(n, (basisSize - i)*n)];
end
[e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18] = deal(e{:});

%--------------------------------------------------------------------------
% Delay-independent composite vectors
%--------------------------------------------------------------------------
es=(B*e1+Bd*e3);

eg=[e2;e3;e4;e6;e7;e12;e13;e15;e16;e17;e18];

c01=[e2;e12;e3;e13];
c02=[e3;e13;e4;e14];
c03=[e1;es-e1;e2;e12];

c1=[e1; es-e1];
c2=[e2; e12];
c3=[e3; e13];
c4=[e4; e14];
c5=[e11;e1-e11;es-2*e1+e11];

E01=[e6 ; e0 ; e0 ; 2*e9-e6 ; e0 ; e0];
E02=[e7 ; e0 ; e0 ; 2*e10-e7 ; e0 ; e0];

E11=[e6-e2 ; e2-e3 ; e12-e13 ; 4*e9-e6-e2 ; e2+e3-2*e6 ; e12+e13-2*e15];
E12=[e7-e3 ; e3-e4 ; e13-e14 ; 4*e10-e7-e3 ; e3+e4-2*e7 ; e13+e14-2*e16];

F11=[e2-e3 ; e12-e13 ; e2+e3-2*e6 ; e12+e13-2*e15 ; e2-e3+6*e6-12*e9 ; e12-e13+6*e15-12*e17];
F12=[e3-e4 ; e13-e14 ; e3+e4-2*e7 ; e13+e14-2*e16 ; e3-e4+6*e7-12*e10 ; e13-e14+6*e16-12*e18];
F1=[F11;F12];

sigma=[-1/2*eye(n) zeros(n) -1/2*eye(n) zeros(n) zeros(n) zeros(n);
    eye(n) zeros(n) zeros(n) zeros(n) zeros(n) zeros(n);
    zeros(n) eye(n) zeros(n) zeros(n) zeros(n) zeros(n);
    -1/6*eye(n) zeros(n) -1/2*eye(n) zeros(n) -1/3*eye(n) zeros(n);
    zeros(n) zeros(n) eye(n) zeros(n) zeros(n) zeros(n);
    zeros(n) zeros(n) zeros(n) eye(n) zeros(n) zeros(n)];

d1=15;    % Lower delay bound d_l   % Set d1=15,25,35,45,50 for example1 and d1=2,6,10,15,20,25 for example2
d2=d1;    % Upper delay bound d_u   % You can select a different dataset for d2 to reduce processing time (with d2 >= d1)
%--------------------------------------------------------------------------
% MAUB Search Loop
%-------------------------------------------------------------------------
primal=1;
while(min(primal)>0)
    d2=d2+1;    % initialise; first action in loop is d2 = d2+1
    d21=d2-d1;

    q1=d21+1;  % q_1(d_21)
    q2=d21+2;  % q_2(d_21)

    %----------------------------------------------------------------------
    % Delay-dependent composite vectors (at d(k)=d_l=d1 and d(k)=d_u=d2)
    %----------------------------------------------------------------------
    % Upsilon_1 : extracts eta_11(k+1)
    Ups1_d1=[e1 ; es ; e2+e12 ; e3+e13 ; e4+e14 ; (d1+1)*e5-e2 ; e6+q1*e7-e3-e4 ; (d1+1)*(d1+2)*e8-(d1+1)*e5 ; 2*e9+q1*q2*e10+(d21-1)*e6-q1*e7-d21*e3];
    Ups1_d2=[e1 ; es ; e2+e12 ; e3+e13 ; e4+e14 ; (d1+1)*e5-e2 ; q1*e6+e7-e3-e4 ; (d1+1)*(d1+2)*e8-(d1+1)*e5 ; q1*q2*e9+2*e10-q1*e6-e7];

    % Gamma_1 : extracts eta_11(k)
    Gam1_d1=[e11 ; e1 ; e2 ; e3 ; e4 ; (d1+1)*e5-e1 ; e6+q1*e7-e2-e3 ; (d1+1)*(d1+2)*e8-(d1+1)*e1 ; 2*e9+q1*q2*e10+d21*e6-q1*e2-q1*e3];
    Gam1_d2=[e11 ; e1 ; e2 ; e3 ; e4 ; (d1+1)*e5-e1 ; q1*e6+e7-e2-e3 ; (d1+1)*(d1+2)*e8-(d1+1)*e1 ; q1*q2*e9+2*e10-q1*e2-e3];

     % Upsilon_2 : extracts eta_12(k+1)
    Ups2_d1=[es-e1 ; es-e2-e12 ; e15+q1*e16-e13-e14 ; d1*es+e2-(d1+1)*e5 ; 2*e17+q1*q2*e18+(d21-1)*e15-q1*e16-d21*e13];
    Ups2_d2=[es-e1 ; es-e2-e12 ; q1*e15+e16-e13-e14 ; d1*es+e2-(d1+1)*e5 ; q1*q2*e17+2*e18-q1*e15-e16];

    % Gamma_2 : extracts eta_12(k)
    Gam2_d1=[e1-e11 ; e1-e2 ; e15+q1*e16-e12-e13 ; (d1+1)*(e1-e5) ; 2*e17+q1*q2*e18+d21*e15-q1*e12-q1*e13];
    Gam2_d2=[e1-e11 ; e1-e2 ; q1*e15+e16-e12-e13 ; (d1+1)*(e1-e5) ; q1*q2*e17+2*e18-q1*e12-e13];

    % F2, F3 (depend on d1 only — recomputed each iter for clarity)
    F2=[(d1+1)*e5-e1 ; e1-e2 ; es-e1-e12 ; 2*(d1+2)*e8-(d1+1)*e5-e1 ; e1+e2-2*e5 ; (d1-1)/(d1+1)*es-e1+e12+2/(d1+1)*e2];

    F3=[e1-e2 ; es-e1-e12 ; e1+e2-2*e5 ; (d1-1)/(d1+1)*es-e1+e12+2/(d1+1)*e2 ; e1-e2+6*e5-12*e8...
        ; ((d1-1)*(d1-2))/((d1+1)*(d1+2))*es-e1-e12-6/(d1+1)*e2+12/(d1+2)*e5];

    %----------------------------------------------------------------------
    % Decision Variables
    %----------------------------------------------------------------------
    P1 = sdpvar(9*n,9*n,'symmetric');
    P2 = sdpvar(5*n,5*n,'symmetric');
    Q1 = sdpvar(2*n,2*n,'symmetric');
    Q2 = sdpvar(2*n,2*n,'symmetric');
    Q3 = sdpvar(2*n,2*n,'symmetric');
    X1 = sdpvar(6*n,6*n,'full');
    X2 = sdpvar(6*n,6*n,'full');
    G1 = sdpvar(11*n,n);
    G2 = sdpvar(11*n,n);
    G3 = sdpvar(11*n,n);
    G4 = sdpvar(11*n,n);

    T1 = sdpvar(n,n,'symmetric');
    T2 = sdpvar(n,n,'symmetric');
    T3 = sdpvar(n,n,'symmetric');
    Z1 = sdpvar(n,n,'symmetric');
    Z2 = sdpvar(n,n,'symmetric');
    Z3 = sdpvar(n,n,'symmetric');

    R11=sdpvar(n,n,'symmetric');
    R22=sdpvar(n,n,'symmetric');
    R33=sdpvar(n,n,'symmetric');
    R12=sdpvar(n,n,'full');
    R13=sdpvar(n,n,'full');
    R23=sdpvar(n,n,'full');

    S11=sdpvar(n,n,'symmetric');
    S22=sdpvar(n,n,'symmetric');
    S33=sdpvar(n,n,'symmetric');
    S12=sdpvar(n,n,'full');
    S13=sdpvar(n,n,'full');
    S23=sdpvar(n,n,'full');

    % \tilde{Z}_i = blkdiag(Z_i,Z_i) for i=1,2,3
    Z1t=blkdiag(Z1,Z1);
    Z2t=blkdiag(Z2,Z2);
    Z3t=blkdiag(Z3,Z3);

    % \Pi_2(\tilde{Z}_i) = blkdiag(\tilde{Z}_i, 3\tilde{Z}_i, 5\tilde{Z}_i) 
    Z1_pi2=blkdiag(Z1t,3*Z1t,5*Z1t);
    Z2_pi2=blkdiag(Z2t,3*Z2t,5*Z2t);
    Z3_pi2=blkdiag(Z3t,3*Z3t,5*Z3t);


    % R_ia and S_a (matrix-separation decomposition, Lemma 3)
    R=[R11 R12 R13;R12' R22 R23;R13' R23' R33];
    R1a=[R11 R12-T1 R13;R12'-T1' R22-T1-Z1 R23-T1;R13' R23'-T1' R33-T1-Z1];
    R2a=[R11 R12-T2 R13;R12'-T2' R22-T2-Z2 R23-T2;R13' R23'-T2' R33-T2-Z2];

    S=[S11 S12 S13;S12' S22 S23;S13' S23' S33];
    Sa=[S11 S12-T3 S13;S12'-T3' S22-T3-Z3 S23-T3;S13' S23'-T3' S33-T3-Z3];

     % Pi_1(.) = blkdiag(., 3.)   [2 blocks]
    R1a_pi1=blkdiag(R1a,3*R1a);
    R2a_pi1=blkdiag(R2a,3*R2a);
    Sa_pi1=blkdiag(Sa,3*Sa);

    % T_i matrices (cal_T_i = blkdiag(T_i, T_i, -T_i, -T_i)) for i=1,2,3
    t1 = blkdiag( T1,  T1, -T1, -T1);
    t2 = blkdiag( T2,  T2, -T2, -T2);
    t3 = blkdiag( T3,  T3, -T3, -T3);
   
    %----------------------------------------------------------------------
    % Psi Terms (Lemma 3 bounds applied to J1, J2, J3)
    %----------------------------------------------------------------------
        
    % Psi_11: at d(k)=d_l -> d_k1=0; at d(k)=d_u -> d_k1=d21
    PSI11_d1=(E01'*R1a_pi1*E11+E11'*R1a_pi1*E01)+c01'*t1*c01;
    PSI11_d2=d21*E01'*R1a_pi1*E01+(E01'*R1a_pi1*E11+E11'*R1a_pi1*E01)+c01'*t1*c01;

    % Psi_12: at d(k)=d_l -> d_k2=d21; at d(k)=d_u -> d_k2=0
    PSI12_d1=d21*E02'*R2a_pi1*E02+(E02'*R2a_pi1*E12+E12'*R2a_pi1*E02)+c02'*t2*c02;
    PSI12_d2=(E02'*R2a_pi1*E12+E12'*R2a_pi1*E02)+c02'*t2*c02;

    PSI21=sigma'*R1a_pi1*sigma+Z1_pi2;
    PSI22=sigma'*R2a_pi1*sigma+Z2_pi2;

    PSI3=1/d1*(F2'*Sa_pi1*F2+F3'*Z3_pi2*F3)+c03'*t3*c03;

    % Psi_4 via Lemma 2 (reciprocally convex bound), alpha=d_k1/d21
    %   alpha=0 at d_l:
    %   alpha=1 at d_u:
    PSI2_d1=[2*PSI21 X1;X1' PSI22];   %alpha = 0
    PSI2_d2=[PSI21 X2;X2' 2*PSI22];   %alpha = 1
    PSI4_d1=1/d21*F1'*PSI2_d1*F1;    % F1 = [F11;F12];
    PSI4_d2=1/d21*F1'*PSI2_d2*F1;    % F1 = [F11;F12];

    PSI_d1=PSI11_d1+PSI12_d1+PSI3+PSI4_d1;
    PSI_d2=PSI11_d2+PSI12_d2+PSI3+PSI4_d2;

    %----------------------------------------------------------------------
    % Phi Terms (LKF forward difference components)
    %----------------------------------------------------------------------
    PHI1_d1=Ups1_d1'*P1*Ups1_d1+Ups2_d1'*P2*Ups2_d1-(Gam1_d1'*P1*Gam1_d1+Gam2_d1'*P2*Gam2_d1);
    PHI1_d2=Ups1_d2'*P1*Ups1_d2+Ups2_d2'*P2*Ups2_d2-(Gam1_d2'*P1*Gam1_d2+Gam2_d2'*P2*Gam2_d2);


    PHI2=c1'*(Q1+q1*Q3)*c1+c2'*(Q2-Q1)*c2-c3'*Q3*c3-c4'*Q2*c4;

    PHI3=c5'*(d21*R+d1*S)*c5;

    % Zero-value terms (Phi_41 to Phi_44) at d(k)=d_l
    z1_d1 = eg'*G1*(e12+e2-e3-e15);         % d_k1=0 -> q_1(d_k1)=1
    z2_d1 = eg'*G2*(e13+e3-e4-q1*e16);      % d_k2=d21 -> q_1(d_k2)=q1
    z3_d1 = eg'*G3*(e12+e2-e6-2*e17);       % d_k1=0 -> q_2(d_k1)=2
    z4_d1 = eg'*G4*(e13+e3-e7-q2*e18);      % d_k2=d21 -> q_2(d_k2)=q2

    % Zero-value terms at d(k)=d_u
    z1_d2 = eg'*G1*(e12+e2-e3-q1*e15);      % d_k1=d21 -> q_1(d_k1)=q1
    z2_d2 = eg'*G2*(e13+e3-e4-e16);         % d_k2=0 -> q_1(d_k2)=1
    z3_d2 = eg'*G3*(e12+e2-e6-q2*e17);      % d_k1=d21 -> q_2(d_k1)=q2
    z4_d2 = eg'*G4*(e13+e3-e7-2*e18);       % d_k2=0 -> q_2(d_k2)=2

    PHI4_d1=z1_d1+z1_d1'+z2_d1+z2_d1'+z3_d1+z3_d1'+z4_d1+z4_d1';
    PHI4_d2=z1_d2+z1_d2'+z2_d2+z2_d2'+z3_d2+z3_d2'+z4_d2+z4_d2';

    PHI_d1=PHI1_d1+PHI2+PHI3+PHI4_d1;
    PHI_d2=PHI1_d2+PHI2+PHI3+PHI4_d2;


    %----------------------------------------------------------------------
    % Full Xi LMI Blocks  (Theorem 1, eq. 15)
    %----------------------------------------------------------------------

    XI_d1=[PHI_d1-PSI_d1 F11'*X2 ; X2'*F11 -d21*PSI22];
    
    XI_d2=[PHI_d2-PSI_d2 F12'*X1' ; X1*F12 -d21*PSI21];


    %----------------------------------------------------------------------
    % LMI Constraints
    %----------------------------------------------------------------------

    LMIs = [P1>=0, P2>=0, Q1>=0, Q2>=0, Q3>=0, R>=0, S>=0, Sa>=0, Z1>=0, Z2>=0, Z3>=0, R1a>=0, R2a>=0, ...
    PSI21>=0, PSI22>=0, XI_d1<=0, XI_d2<=0];

    %----------------------------------------------------------------------
    % Solve
    %----------------------------------------------------------------------
    
    option = sdpsettings('solver', 'sdpt3', 'verbose', 0);
    sol = optimize(LMIs, [], option);

    [primal, ~] = check(LMIs);

    if ((sol.problem == 0))
        if min(primal) <= 0
            disp('Infeasible problem');
        else
            disp('Successfully solved');
            disp(['Delay upper bound d_u = ', mat2str(d2)]);
            % check(LMIs)
        end
    else
        yalmiperror(sol.problem);
        break;
    end
end          % <-- while loop end
%--------------------------------------------------------------------------
% Report MAUB
%--------------------------------------------------------------------------
disp(['Delay upper bound d_u = ', mat2str(d2 - 1)]);
