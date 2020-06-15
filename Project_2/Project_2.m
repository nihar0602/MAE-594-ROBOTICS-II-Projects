clear all; clc; close all;
X=10;
Y=0;
Theta_2=atan(Y/X);
R0_6=[1 0 0; 0 1 0; 0 0 1];
R0_3=[-sin(Theta_2) 0 cos(Theta_2); 
    cos(Theta_2) 0 sin(Theta_2);
    0 1 0];
R3_6 = inv(R0_3)*R0_6
Theta_5 = acos(R3_6(3,3))
Theta_6 = acos(R3_6(3,1)/sin(Theta_5))
Theta_4 = acos(R3_6(2,3)/sin(Theta_5))