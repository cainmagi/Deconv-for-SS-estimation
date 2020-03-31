function [data] = script_cut_data(s, subjects, starttime, endtime)
%script_cut_data
%   Cut the required parts for the following processing.
%   This step aims at reducing the calculation cost.
%   Argument:
%       s: the dataset s.
%       subjects: a list of indices, which parts of the data would be
%           cut for the following processing.
%       starttime: the start time of cut window.
%       endtime: the end time of cut window.
%   Returns:
%       data: the cut data part, which would be used further.
if nargin < 2
    subjects = 11;
end
if nargin < 3
    starttime = 200;
end
if nargin < 4
    endtime = starttime + 200;
end
for sub = subjects
    Fsy = 1;
    Fsu = 4;
    Fsq = 5;
    Tsy = 1/Fsy;

    ub = [1.4 6.0 1.0]';
    lb = [0.1 1.5 0.01]';
    
    y1_ton_phas = s(sub).sc(2).yproc;
    y2_ton_phas = s(sub).sc(3).yproc;
    y3_ton_phas = s(sub).sc(1).yproc;
    
    Y = [y1_ton_phas y2_ton_phas y3_ton_phas];
    
    Fss = 100;
    [sigma] = get_sigma_nchannel(Y((starttime*Fss):(endtime*Fss),:), Fss);
    
    [~, c] = size(Y);
    lag = [];
    lag(1) = 0; % first lag is zero, here the reference signal is the channel one.
    
    for i = 2:c
        y_ref_temp = Y(:,1);
        y_com_temp = Y(:,i);

        C = xcorr(y_ref_temp,y_com_temp);

        [~, idx] = max(C);

        lag(i) = round(length(C)/2-idx+1);
        y_comp = [y_com_temp(max(lag(i)-1,1):end); y_com_temp(end)*ones(max(lag(i)-2,0),1)];
        Y(:,i) = y_comp;
    end
    
        
    % per channel filtering and downsampling the phasic component
    Y_ = [];
    for i = 1:c
       Y(:,i) = LowPassFilter(Y(:,i), Fss, 0.5, 64);
       Y_(:,i) = downsample(Y(:,i), Fss/Fsy);
    end
    Y = Y_;
    
    [Ny,~] = size(Y);
    ty = 0:Tsy:(Ny-1)*Tsy;

    % Take a window of the signal
    lls = starttime; % window start time in second
    rrs = endtime; % window end time in second
    ll = lls*Fsy;
    rr = rrs*Fsy;
    
    ty = ty(ll:rr)';
    
    Y = Y(ll:rr,:);

    % make the data structure for the function
    
    data.y = Y;
    data.sigma = sigma;
    data.ub = ub;
    data.lb = lb;
    data.Fsu = Fsu;
    data.Fsy = Fsy;
    data.Fsq = Fsq;
    data.minimum_peak_distance = 20;
   
end
end

