%preprocessing
%---------------------------------------------------
% artifacts removal with heuristic approach
% denoising, tonic part removal
%---------------------------------------------------

clear all;
clc;
close all;

load('auditory_stimulation_dataset\s.mat');

for i = 1:26
    chMax = 3;
    tic
    figure,
    for ch = 1:chMax
            subplot(chMax,1,ch);
            y = s(i).sc(ch).y;
            Fss = s(i).sc(ch).sampling_rate;
            N = length(y);
            ty = (0:length(y)-1)/Fss;
            plot(ty, y);
            [~, idxs1]=findpeaks(diff(y),'MinPeakProminence',0.1);
            [~, idxs2]=findpeaks(-diff(y),'MinPeakProminence',0.1);
            hold on;
            plot(ty(idxs1),y(idxs1),'ro');
            plot(ty(idxs1),y(idxs1),'go');
    %%      noise spikes removal using spline interpolation

            y = [y(75)*ones(75,1); y(76:end-75); y(end-74)*ones(75,1)];

            for k = 1:length(idxs1)
                if(idxs1(k) >100 && idxs1(k) < length(y)-100)
                    x = [idxs1(k)-50 idxs1(k)-25 idxs1(k)+25 idxs1(k)+50];
                    v = y(x);
                    xq = [idxs1(k)-50:idxs1(k)+50];
                    vq = interp1(x,v,xq,'spline');
                    y(xq) = vq;
                end
            end
            for k = 1:length(idxs2)
                if(idxs2(k) >300 && idxs2(k) < length(y)-300)
                    x = [idxs2(k)-50 idxs2(k)-25 idxs2(k)+25 idxs2(k)+50];
                    v = y(x);
                    xq = [idxs2(k)-50:idxs2(k)+50];
                    vq = interp1(x,v,xq,'spline');
                    y(xq) = vq;
                end
            end

            %% LowPassFilter
            y = LowPassFilter(y, Fss, 3, 64);
            hold on,plot(ty,y,'r');

            [r, p, tt, l, d, e, obj] = cvxEDA(y, 1/Fss, 2, 0.7, 6, 8e-4, 1e-2, 'quadprog');
            hold on,stem(ty,p/10,'g.');
            s(i).sc(ch).phasic = r;
            s(i).sc(ch).tonic = tt;
            s(i).sc(ch).driver = p;
            ss(i) = s(i);
            hold on,plot(ty,tt,'r--');
            hold on,plot(ty,r,'b--');
    end
    savefig(['tonic_phasic_subject_',num2str(i),'.fig']);
    close all;
    toc
end

save('ss_8_14_18.mat','ss');