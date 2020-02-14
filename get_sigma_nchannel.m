function [sigma] = get_sigma_nchannel(Y, Fss)

    [r,c] = size(Y);

    for i=1:c
        y = Y(:,i);
        yL = LowPassFilter(y, Fss, 3, 64);
        yH = y-yL;
        sigma(i) = sqrt(var(yH)*(Fss/2)/(Fss/2 - 3));   
    end

end
