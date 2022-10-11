clear; close all; clc;
h_ba = 1.05457266e-34;
gama = 2*pi*42.58e6;
T = 310;
b0 = 3;
kb = 1.38064852e-23;
delta_e = h_ba*gama*b0;
result = exp(-delta_e/(kb*T));
ratio = (1-result)/(1+result);