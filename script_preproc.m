function [s] = script_preproc(s)
%Data preprocessing.
%   Input the dataset (struct), and returns the pre-processed results.
%   Argument:
%       s: the dataset s.
%   Returns:
%       s: the dataset s, where the y vectors are added with processed
%       ones.
chMax = numel(s(1).sc); % 3
for i = 1:numel(s) % 26
    for ch = 1:chMax
        y = s(i).sc(ch).y;
        Fss = s(i).sc(ch).sampling_rate;
        [~, idxs1]=findpeaks(diff(y),'MinPeakProminence',0.1);
        [~, idxs2]=findpeaks(-diff(y),'MinPeakProminence',0.1);
        % noise spikes removal using spline interpolation
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
        
        % LowPassFilter
        y = LowPassFilter(y, Fss, 3, 64);
        s(i).sc(ch).yproc = y;
    end
end
end