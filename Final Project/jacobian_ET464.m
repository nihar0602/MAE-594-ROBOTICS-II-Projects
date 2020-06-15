function [F, J] = jacobian_ET464(Q,x,y,z)
theta1 = Q(1);
theta2 = Q(2);


DH02 = [ cos(theta1)*cos(theta2), -cos(theta1)*sin(theta2),  sin(theta1), (13*cos(theta1)*cos(theta2))/50; cos(theta2)*sin(theta1), -sin(theta1)*sin(theta2), -cos(theta1), (13*cos(theta2)*sin(theta1))/50 ;             sin(theta2),              cos(theta2),            0,    (13*sin(theta2))/50 + 33/100 ;                       0,                        0,            0,                               1];
DH01 = [ cos(theta1), 0,  sin(theta1),      0; sin(theta1), 0, -cos(theta1),      0;            0, 1,            0, 33/100;           0, 0,            0,      1];
% q3 = Q(3);



F = [DH02(1,4)-x ; ...
    DH02(2,4)-y; ...
    DH02(3,4)-z];
J = [cross([0 0 1]', DH02(1:3,4)-[0 0 0]') cross(DH01(1:3,3), DH02(1:3,4)-DH01(1:3,4));
            [0 0 1]' DH01(1:3,3)];

% J = [-sin(q1)-sin(q1+q2)-sin(q1+q2+q3) , -sin(q1+q2)-sin(q1+q2+q3), -sin(q1+q2+q3);...
%     cos(q1)+cos(q1+q2)+cos(q1+q2+q3) , cos(q1+q2)+cos(q1+q2+q3),cos(q1+q2+q3)];