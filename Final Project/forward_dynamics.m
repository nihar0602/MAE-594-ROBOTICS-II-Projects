clear all;
close all;
clc;

syms theta1 theta2 %theta3(t) theta4(t) theta5(t) theta6(t) %(t) theta2(t) theta3(t) theta4(t) theta5(t) theta6(t)%[1 6] % symbolic variables
syms alpha r d % symbolic variables
syms qdot1 qdot2 %qdot3(t) qdot4(t) qdot5(t) qdot6(t)%(theta1,t) qdot2(theta2,t) qdot3(theta3,t) qdot4(theta4,t) qdot5(theta5,t) qdot6(theta6,t)
syms qddot1 qddot2 %qddot3(t) qddot4(t) qddot5(t) qddot6(t) 
DH = [cos(theta1) -sin(theta1)*cos(alpha) sin(theta1)*sin(alpha) r*cos(theta1); % generalized DH form
    sin(theta1) cos(theta1)*cos(alpha) -cos(theta1)*sin(alpha) r*sin(theta1);
    0 sin(alpha) cos(alpha) d;
    0 0 0 1];
DH01 = subs(DH, [alpha r d], [pi/2 0 0.330]); % DH b/w frame 0 and 1
DH12 = subs(DH, [theta1 alpha r d], [theta2 0 0.260 0]); % DH b/w frame 1 and 2
% DH23 = subs(DH, [theta1 alpha r d], [theta3 + pi/2 pi/2 0 0]); % DH b/w frame 2 and 3
% DH34 = subs(DH, [theta1 alpha r d], [theta4 + pi pi/2 0 290]); % DH b/w frame 3 and 4
% DH45 = subs(DH, [theta1 alpha r d], [theta5 + pi pi/2 0 0]); % DH b/w frame 4 and 5
% DH56 = subs(DH, [theta1 alpha r d], [theta6 0 0 70]); % DH b/w frame 5 and 6
DH02 = DH01*DH12; % DH b/w frame 0 and 2
% DH03 = DH02*DH23; % DH b/w frame 0 and 3
% DH04 = DH03*DH34; % DH b/w frame 0 and 4
% DH05 = DH04*DH45; % DH b/w frame 0 and 5
P01 = subs(DH, [alpha r d], [pi/2 0 0.165]); % Substituting DH values in case of cg at 1
P12 = subs(DH, [theta1 alpha r d], [theta2 0 0.130 0]); % Substituting DH values in case of cg at 1
% P23 = subs(DH, [theta1 alpha r d], [theta3 + pi/2 pi/2 0 0]); % Substituting DH values in case of cg at 1
% P34 = subs(DH, [theta1 alpha r d], [theta4 + pi pi/2 0 145]); % Substituting DH values in case of cg at 1
% P45 = subs(DH, [theta1 alpha r d], [theta5 + pi pi/2 0 0]); % Substituting DH values in case of cg at 1
% P56 = subs(DH, [theta1 alpha r d], [theta6 0 0 35]); % Substituting DH values in case of cg at 1
CG1 = P01;
CG2 = DH01*P12;
% CG3 = DH01*DH12*DH23*P34;
% CG4 = DH01*DH12*DH23*DH34*DH45*P56;
JACCG1 = [cross([0 0 1]', CG1(1:3,4)-[0 0 0]');
            [0 0 1]'];
JACCG2 = [cross([0 0 1]', CG2(1:3,4)-[0 0 0]') cross(DH01(1:3,3), CG2(1:3,4)-DH01(1:3,4));
            [0 0 1]' DH01(1:3,3)];
% JACCG3 = [cross([0 0 1]', CG3(1:3,4)-[0 0 0]') cross(DH01(1:3,3), CG3(1:3,4)-DH01(1:3,4)) cross(DH02(1:3,3), CG3(1:3,4)-DH02(1:3,4)) cross(DH03(1:3,3), CG3(1:3,4)-DH03(1:3,4));
%             [0 0 1]' DH01(1:3,3) DH02(1:3,3) DH03(1:3,3)];            
% JACCG4 = [cross([0 0 1]', CG4(1:3,4)-[0 0 0]') cross(DH01(1:3,3), CG4(1:3,4)-DH01(1:3,4)) cross(DH02(1:3,3), CG4(1:3,4)-DH02(1:3,4)) cross(DH03(1:3,3), CG4(1:3,4)-DH03(1:3,4)) cross(DH04(1:3,3), CG4(1:3,4)-DH04(1:3,4)) cross(DH05(1:3,3), CG4(1:3,4)-DH05(1:3,4));
            %[0 0 1]' DH01(1:3,3) DH02(1:3,3) DH03(1:3,3) DH04(1:3,3) DH05(1:3,3)]; 
m1 = 2710e-9*11913000;
m2 = 2710e-9*9386000;
% m3 = 2710e-9*6097250;
% m4 = 2710e-9*1471750;
I1 = 1/12*m1*[0.190.^2+0.330.^2 0 0;
                            0 0.190.^2+0.330.^2 0;
                            0 0 0.190.^2+0.190.^2];
I2 = 1/12*m2*[0.190^2+0.190^2 0 0; 0 0.190^2+0.260^2 0; 0 0 0.190^2+0.260^2];
% I3 = 1/12*m3*[145.^2+290.^2 0 0;
%                            0 145.^2+290.^2 0;
%                            0 0 145.^2+145.^2];
% I4 = 1/12*m4*[145.^2+70.^2 0 0;
%                            0 145.^2+70.^2 0;
%                            0 0 145.^2+145.^2];

omega1 = JACCG1(4:6,:)*qdot1;
omega2 = DH01(1:3,1:3)*JACCG2(4:6,:)*[qdot1; qdot2];
% omega3 = DH03(1:3,1:3)*JACCG3(4:6,:)*[qdot1; qdot2; qdot3; qdot4];
% omega4 = DH05(1:3,1:3)*JACCG4(4:6,:)*[qdot1; qdot2; qdot3; qdot4; qdot5; qdot6];
%% Kinetic Energy 
D1 = m1*(JACCG1(1:3,:))'*(JACCG1(1:3,:))*qdot1 + (JACCG1(4:6,:))'*I1*omega1;
K_E1 = 0.5 * qdot1 * D1
D2 = m2*(JACCG2(1:3,:))'*(JACCG2(1:3,:))*[qdot1;qdot2] + (JACCG2(4:6,:))'*DH01(1:3,1:3)*I2*omega2;
K_E2 = 0.5 * [qdot1; qdot2]' * D2;
% D3 = m3*(JACCG3(1:3,:))'*(JACCG3(1:3,:))*[qdot1;qdot2; qdot3; qdot4] + (JACCG3(4:6,:))'*DH03(1:3,1:3)*I3*omega3;
% K_E3 = 0.5 * [qdot1;qdot2; qdot3; qdot4]' * D3;
% D4 = m4*(JACCG4(1:3,:))'*(JACCG4(1:3,:))*[qdot1;qdot2; qdot3; qdot4; qdot5; qdot6] + (JACCG4(4:6,:))'*DH05(1:3,1:3)*I4*omega4;
% K_E4 = 0.5 * [qdot1;qdot2; qdot3; qdot4; qdot5; qdot6]' * D4 % 
%% Potential Energy
P_E1 = [0 0 9.800]*CG1(1:3,4)*m1; % Potential Enegry for link 1
P_E2 = [0 0 9.800]*CG2(1:3,4)*m2; % Potential Enegry for link 2
% P_E3 = [0 0 9800]*CG3(1:3,4)*m3; % Potential Enegry for link 3
% P_E4 = [0 0 9800]*CG4(1:3,4)*m4; % Potential Enegry for link 4
%% Lagrangian
%L = (K_E1 + K_E2) - (P_E1 + P_E2) % Lagrangian
K_E = K_E1 + K_E2;
P_E = P_E1 + P_E2;
KEwrtqdot1 = diff(K_E,qdot1);
KEwrtqdot2 = diff(K_E,qdot2);
KEwrttheta1 = diff(K_E,theta1);
KEwrttheta2 = diff(K_E,theta2);
KEwrtqtt1 = diff(KEwrtqdot1,theta1)*qdot1 + diff(KEwrttheta1,qdot1)*qddot1 + diff(KEwrtqdot1,theta2)*qdot2 + diff(KEwrttheta1,qdot2)*qddot2
KEwrtqtt2 = diff(KEwrtqdot2,theta1)*qdot1 + diff(KEwrttheta2,qdot1)*qddot1 + diff(KEwrtqdot2,theta2)*qdot2 + diff(KEwrttheta2,qdot2)*qddot2
PEwrttheta1 = diff(P_E,theta1);
PEwrttheta2 = diff(P_E,theta2);

%%
% Tau1 = 700000000;

% Tau2 = 600000000;
syms Tau1 Tau2
eom1 = simplify(Tau1 - (KEwrtqtt1 - KEwrttheta1 + PEwrttheta1))
eom2 = simplify(Tau2 - (KEwrtqtt2 - KEwrttheta2 + PEwrttheta2))

Sol = solve(eom1,eom2,qddot1,qddot2);
Sol.qddot1 = simplify(Sol.qddot1);
Sol.qddot2 = simplify(Sol.qddot2);


syms y1 y2 y3 y4 y5 y6
fx1=subs(Sol.qddot1,[theta1,qdot1,theta2,qdot2],[y1,y2,y3,y4])
fx2=subs(Sol.qddot2,[theta1,qdot1,theta2,qdot2],[y1,y2,y3,y4])

% diary ('fx1.txt');
% fx1 % your symbolic variable
% diary off
% 
% diary ('fx2.txt');
% fx2 % your symbolic variable
% diary off


%%
close all; clear; clc;
load ws.mat
init_t=0;
final_t=1;
dt=0.001;
N= (final_t-init_t)/dt;
t_span=linspace(init_t,final_t,N);
% load theta.mat
% t_span=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]

x0=[0.1 0.1 pi/3 0.1]';
% u=[3 0]';

x=zeros(4,N);
x(:,1)=x0;

% for i=1:N-1
%     x(:,i+1) = two_dof_arm_dyn(x(:,i),u)*dt+x(:,i);
% end

% [t,y] = ode45(@two_dof_arm_dyn_for_ODE,[init_t final_t],x0);
% [t,y] = ode45(@two_dof_arm_dyn_for_ODE,t_span,x0);
Tau = [1000;1000];
for i=2:N
    x(:,i) = RK4(@DE, dt, x(:,(i-1)),Tau);
end

% 
% L1=1;
% L2=1;
% M1=1;
% M2=1;
% g=9.8;


x = x'

mov_cnt = 1;
figure; hold on;
for i=1:N-1
%     if(mod(i,50)==1)
        
        clf;
        x1=x(i,1);
        x2=x(i,3);
        p0x=0;
        p0y=0;
        p0z=0;
        PH01(1,4) = subs(DH01(1,4), [theta1, theta2],[x1, x2]);
        PH01(2,4) = subs(DH01(2,4), [theta1, theta2],[x1, x2]);
        PH01(3,4) = subs(DH01(3,4), [theta1, theta2],[x1, x2]);
        PH02(1,4) = subs(DH02(1,4), [theta1, theta2],[x1, x2]);
        PH02(2,4) = subs(DH02(2,4), [theta1, theta2],[x1, x2]);
        PH02(3,4) = subs(DH02(3,4), [theta1, theta2],[x1, x2]);
        p1x = PH01(1,4);
        p1y = PH01(2,4);
        p1z = PH01(3,4);
        p2x = PH02(1,4);
        p2y = PH02(2,4);
        p2z = PH02(3,4);
        px=[p0x p1x p2x];
        py=[p0y p1y p2y];
        pz=[p0z p1z p2z];
        disp(px)
        plot3(px,py,pz,'ro-');
        axis([-1 1 -1 1 -1 1]);
        pause(0.001);
        MM(mov_cnt)=getframe;
        mov_cnt=mov_cnt+1;
    
end

movie(MM)
movie2avi(MM,'twoDOFarm.avi');
%%

