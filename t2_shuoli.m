clear; close all; clc;

m_0 = 1;
t_1 = 500; 
t_2 = 50;
angle_flip = 10*pi/180; % flip angle 
t_r = 50; % repetition time 
pulse_num = 50; % number of pulses
pulse_length = 100; % the length of a single pulse
M_z = zeros(pulse_num);
M_xy = zeros(pulse_num);
M_z_2 = zeros(pulse_length,pulse_num);
M_xy_2 = zeros(pulse_length,pulse_num);


M_z_before_pulse = m_0; % initialization of Mz
M_xy_before_pulse = 0; % initialization of Mxy
t = 0:1:(pulse_length-1);
t = t.*t_r./(pulse_length-1);
for p = 1:1:pulse_num
    M_z(p) = M_z_before_pulse;
    Mz_after_pulse = M_z_before_pulse*cos(angle_flip) - M_xy_before_pulse*sin(angle_flip);
    Mxy_after_pulse = M_z_before_pulse*sin(angle_flip) + M_xy_before_pulse*cos(angle_flip);
    M_xy(p) = Mxy_after_pulse;
    M_z_2(:,p) = (Mz_after_pulse - m_0) * exp(-t./t_1) + m_0;
    M_xy_2(:,p) = Mxy_after_pulse * exp(-t./t_2);
    M_z_before_pulse = M_z_2(pulse_length,p);
    M_xy_before_pulse = M_xy_2(pulse_length,p);
end

figure(1);
subplot(2,2,1);
plot(M_z); title('M_z'); 
xlim([0,pulse_num]); ylim([-1,1]);
subplot(2,2,3);
plot(M_xy); title('M_x_y'); 
xlim([0,pulse_num]); ylim([-1,1]);

subplot(2,2,2);
plot(reshape(M_z_2,1,pulse_num*pulse_length)); 
title('M_z'); 
xlim([0,pulse_num*pulse_length]); ylim([-1,1]);
subplot(2,2,4);
plot(reshape(M_xy_2,1,pulse_num*pulse_length)); 
title('M_x_y'); 
xlim([0,pulse_num*pulse_length]); ylim([-1,1]);
