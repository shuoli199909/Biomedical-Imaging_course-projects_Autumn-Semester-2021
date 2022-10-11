clear; close all; clc;

%matrix = 256;
%[filter,y,filter_ori] = CalcFilter(matrix+1,0);
%figure(1);
%plot(filter,'b-','LineWidth',1); hold on;
%plot(y,'r-','LineWidth',1); hold on;
%plot(filter_ori,'b-','LineWidth',1);
%xlim([0 257]);
%ylim([0 1]);
%title('Ideal high-pass filter');
%legend('High-pass filter(ideal)')

%figure(2);
%plot(fftshift(real(fft(fftshift(filter)))),'b-','LineWidth',1); hold on;
%plot(fftshift(real(fft(fftshift(y)))),'r-','LineWidth',1); hold on;
%plot(fftshift(real(fft(fftshift(filter_ori)))),'b-','LineWidth',1);
%xlim([129-20 129+20]);
%ylim([-60 150]);
%title('Ideal high-pass filter');
%legend('High-pass filter(ideal)')

load('phantom_save.mat');
sbp_fft = fft(phantom.sbp);
fbp_fft = fft(phantom.fbp);
ori_fft = fft(phantom.discrete);

diff_fft = fbp_fft - ori_fft;
figure(1);
DisplayData(diff_fft);
title('fft(phantom.fbp-phantom.discrete)');

function [filter,y,filter_ori] = CalcFilter(matrix,treshold)

    % --------------------------------------------------------------------
    % Define high-pass filter according to |u|
    % --------------------------------------------------------------------
    filter = abs(-fix(matrix/2):+fix(matrix/2));    % filter placeholder
    
    % TASK 2.3 FILL IN HERE
    for i = 0:matrix/2+1
        filter(filter==i) = treshold+(1-treshold)*(i/(matrix/2));
        %filter(filter==i) = 0.5 - 0.5*cos(i*pi/180);
    end
    filter_ori = filter;
    x = 1:length(filter);
    a = 6.05e-5;
    y = 1-a*(x-0.5*length(filter)).^2;
    filter = filter.*y;
    %filter = filter/max(filter);
end
