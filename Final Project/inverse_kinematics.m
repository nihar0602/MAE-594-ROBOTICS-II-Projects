close all;clear all;clc;

%Desired Trajectory
t = 0:0.05:2*pi;
y = 0.260*cos(t);
x = 0.260*sin(t);
z = 0.330 + 0.100*cos(t);

theta1 = zeros(length(t),1);
% q2 = zeros(length(t),1);
theta2 = zeros(length(t),1);
% q3 = zeros(length(t),1);

Q = [-pi/3 ;pi/3];

for i=1:length(t)
    
    while(1)
        [F, J] = jacobian_ET464(Q,x(i),y(i),z(i));
        Qnew = Q-pinv(J(1:3,:))*F;
        
        if norm(Qnew-Q)<0.001
            break
        end
        Q = Qnew;
    end
    theta1(i)=Q(1);
    theta2(i)=Q(2);
    Q =[theta1(i);theta2(i)];
    %Q = [-pi/3 ;pi/3 ;pi/3];

%     XX = [0,cos(q1(i)),cos(q1(i))+cos(q1(i)+q2(i)),cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i))];
%     YY = [0,sin(q1(i)),sin(q1(i))+sin(q1(i)+q2(i)),sin(q1(i))+sin(q1(i)+q2(i))+sin(q1(i)+q2(i)+q3(i))];
    DH01 = [ cos(theta1(i)), 0,  sin(theta1(i)),      0; sin(theta1(i)), 0, -cos(theta1(i)),      0;            0, 1,            0, 33/100;           0, 0,            0,      1];

    DH02 = [ cos(theta1(i))*cos(theta2(i)), -cos(theta1(i))*sin(theta2(i)),  sin(theta1(i)), (13*cos(theta1(i))*cos(theta2(i)))/50; cos(theta2(i))*sin(theta1(i)), -sin(theta1(i))*sin(theta2(i)), -cos(theta1(i)), (13*cos(theta2(i))*sin(theta1(i)))/50 ;             sin(theta2(i)),              cos(theta2(i)),            0,    (13*sin(theta2(i)))/50 + 33/100 ;                       0,                        0,            0,                               1];

    XX = [0,DH01(1,4),DH02(1,4)];
    YY = [0,DH01(2,4),DH02(2,4)];
    ZZ = [0,DH01(3,4),DH02(3,4)];
    cla
    plot3(XX,YY,ZZ,'-o','Markersize', 5, 'Linewidth',1)
    hold on
    plot3(x,y,z,':k')
    xlim([-1 1])
    ylim([-1 1])
    zlim([-1 1])
    pause(0.05)
end





