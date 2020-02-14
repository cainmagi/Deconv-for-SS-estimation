%plot all three channels for asiolomar extention
clear variables;
close all;
clc;

load('ss_8_14_18.mat')

Fss = 100;
for sub =11 %[11,12,14,15,18,19, 20, 21, 23, 25, 26, 7, 8, 10]
    figure('units','normalized','outerposition',[0 0 1 1]),
    subplot(311), plot(ss(sub).sc(1).y(200*Fss:400*Fss));
    subplot(312), plot(ss(sub).sc(2).y(200*Fss:400*Fss));
    subplot(313), plot(ss(sub).sc(3).y(200*Fss:400*Fss));
    
    figure('units','normalized','outerposition',[0 0 1 1]),
    subplot(311), plot(ss(sub).sc(1).phasic(200*Fss:400*Fss));
    subplot(312), plot(ss(sub).sc(2).phasic(200*Fss:400*Fss));
    subplot(313), plot(ss(sub).sc(3).phasic(200*Fss:400*Fss));
end