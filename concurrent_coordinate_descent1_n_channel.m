function [data] = concurrent_coordinate_descent1_n_channel(data)
%structure data has following variables
%number_of_channels, y, sigma, ub, lb, Fsu, Fsy, minimum_peak_distance

%y:                     the matrix containing all the channel data (delay adjusted). Each column corresponds to each channel
%sigma:                 array of all the standard deviation of the noise levels
%ub:                    upper bound for rise time, decay time and the attenuation term 
%lb:                    lower bound for "     "     "
%Fsu:                   desired sampling frequency of the neural stimuli 
%Fsy:                   sampling frequency of the skin conductance observations
%minimum_peak_distance: desired minimum distance between two consecutive impulses

%% approximate number of peaks
    y   = data.y;
    Fsy = data.Fsy;
    Fsu = data.Fsu;
    ub  = data.ub;
    lb  = data.lb;

    [~,locs]=findpeaks(y(:,1), 'MinPeakDistance', data.minimum_peak_distance * data.Fsy);
    ubn=length(locs)+20;

%% random initialization
    tau_j_1(1) = lb(1) + (ub(1)-lb(1))*rand(1,1);   tau_j_1(2) = lb(2) + (ub(2)-lb(2))*rand(1,1);

    data.number_of_channels = min(size(data.y)); %assuming data length is greater than the number of channels

    for i = 1:data.number_of_channels-1
        tau_j_1(2+i) = lb(3) + (ub(3)-lb(3))*rand(1,1);
    end

    Nu = length(y(:,1)) * Fsu/Fsy; 

    Tsy = 1/Fsy; %min
    Tsu = 1/Fsu; %min
    ty = 0:Tsy:(y(:,1)-1)*Tsy;
    tu = 0:Tsu:(Nu-1)*Tsu;

%  define blank arrays
A = [];    B = [];    y_ = [];
%% initialization
for j=1:30
%  fprintf('parallel loop %d: j = %d', parloop,j);
    [A_, B_] = create_A_B_matrix_ss_multires(tau_j_1(1:2), Nu, Fsu, Fsy);
%  define dummy arrays for calculation
    A = [];    B = [];    y_ = [];
    alpha = tau_j_1(3:end);
    for  k = 1:data.number_of_channels
        y_ = [y_; y(:,k) - A_ * [0; y(1,k)]];
        A = [A; A_];
        if( k == 1)
            B = [B; 1*B_];
        else
            B = [B; alpha(k-1)*B_];
        end
    end
%  initialize u with all ones
    uj_1 = ones(Nu,1);

%  run heuristic IRLS algorithm for sparse recovery (CS regime)
    [ uj ] = focuss_modified2(uj_1, y_, B, ubn, true, 15, 0.5, 1e-4, 1e-5 , round(data.minimum_peak_distance * Fsu));
    
    data.u = uj;
    tau_j = SystemID_interior_point1(tau_j_1, data);
    tau_j_1 = tau_j;
end

%% coordinated descent
count = 0;
maxiter = 500;
J = zeros(maxiter,1);
while(1)
    [A_, B_] = create_A_B_matrix_ss_multires(tau_j_1(1:2), Nu, Fsu, Fsy);
%     define dummy arrays for calculation
    A = [];    B = [];    y_ = [];
    alpha = tau_j_1(3:end);
    for  k = 1:data.number_of_channels  
        y_ = [y_; y(:,k) - A_ * [0; y(1,k)]];
        A = [A; A_];
        if( k == 1 )
            B = [B; B_];
        else
            B = [B; alpha(k-1)*B_];
        end
    end
   
    uj_1 = uj;
    
    [uj, J(count+1), Reg] = focussreg3_modified(uj_1, y_, B, 10);
    
    tau_j = SystemID_interior_point1(tau_j_1, data);
    
    %[tau_j,ub,lb]= handle_stagnation_in_boundary(y1,y2, sigma1, sigma2, uj, tau_j, ub, lb,  1e-4, tu, Fsy, Fsu);
    
    tau_j_1 = tau_j;
    alpha_conv_flag = 1;
%     fprintf('count = %d---->   ',count);
%     fprintf('tau1 = %d,   tau2 = %d,   tau3 = %d,   tau4 = %d\n', tau_j(1),tau_j(2),tau_j(3),tau_j(4));
    if(round((uj)*5e1)/5e1 == round((uj_1)*5e1)/5e1)
                disp('convergence achieved for u');
                 if(round(tau_j(1)*1e2)/1e2 == round(tau_j_1(1)*1e2)/1e2 && round(tau_j(2)*1e2)/1e2 == round(tau_j_1(2)*1e2)/1e2 && round(tau_j(2)*1e2)/1e2 == round(tau_j_1(2)*1e2)/1e2)
                     disp('convergence achived for tau1 and tau2');
                     
                     %check convergence for alpha
                     for i = 1:data.number_of_channels-1
                        alpha_conv_flag = (alpha_conv_flag && (round(tau_j(i+2)*1e2)/1e2 == round(tau_j_1(i+2)*1e2)/1e2));
                     end
                     if(alpha_conv_flag == 1)
                        disp('convergence achived for alpha(s) achived');
                        convergenceFlag = 1;
                        break;
                     end
                 end
    end
    count = count+1;
    if(count>maxiter)
        convergenceFlag = 0;
        break;
    end
end


%% reconstruct the waves and save the results in same structure

    tau_j = tau_j(:)';
    [A_, B_] = create_A_B_matrix_ss_multires(tau_j(1:2), Nu, data.Fsu, data.Fsy);
    
    for k = 1:data.number_of_channels
        if(k==1)
            alpha = 1;
        else
            alpha = tau_j(2+k-1);
        end
        y_rec1 = A_*[0;y(1,k)]+alpha*B_*uj;
        y_rec(:,k) = y_rec1;
    end

    data.lambda = Reg;
    data.y_rec = y_rec;
    data.tau_j = tau_j;
    data.uj = uj;
    data.convergenceFlag = convergenceFlag;
    
    data.cost1 = 0;
    data.cost2 = 0;
    data.R_2 = [];
    for i = 1:data.number_of_channels
        data.cost1 = data.cost1 + 0.5*norm(y(:,i)-y_rec(:,i),2)^2/data.sigma(i);
        data.R_2(i) = 1 - var(y(:,i)-y_rec(:,i))/var(y(:,i));
    end
    data.cost2 = data.cost1 + data.lambda * norm(data.uj, 1);
    
end

%% system identification
function tau_est = SystemID_interior_point1(tau_prev, data)
%y1, y2, sigma1, sigma2, ty, u, tu, ub, lb
    lb = data.lb; ub = data.ub;
    lb = lb(:); ub = ub(:); 
    
    for i = 1:data.number_of_channels-2
        lb = [lb; lb(3)]; ub = [ub; ub(3)];
    end
    
    fun = @(x)cost_function_interior_point(x,data);
    %x0 = (ub+lb)/2;  
    x0 = tau_prev(:);
    A = []; b = [];  Aeq = [];  beq = [];
    nonlcon = [];
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off'); % run interior-point algorithm
    tau_est = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    tau_est = tau_est(:)';
end

%% cost function
function J = cost_function_interior_point(tau, data)

    tau = tau(:)';
    Nu = length(data.y(:,1)) * data.Fsu/data.Fsy; 
    [A1, B1] = create_A_B_matrix_ss_multires(tau(1:2), Nu, data.Fsu, data.Fsy);
    
    J = 0;
    for k = 1:data.number_of_channels
        if(k==1)
            alpha = 1;
        else
            alpha = tau(2+k-1);
        end
        J = J + 0.5*norm(data.y(:,k)-A1*[0;data.y(1,k)]-alpha*B1*data.u,2)^2/data.sigma(k)^2;
    end

end


function [A, B] = create_A_B_matrix_ss_multires(tau, Nu, Fsu, Fsy)

% theta, skin
% Nu
% Fsu
% Fsy

    Fs = Fsu;
    Ts = 1/Fs;
    A = [-1/tau(1) 0; +1/tau(2) -1/tau(2)];      B = [1/tau(1) ;0];  C = [0 1];  D = 0;

    sys = ss(A,B,C,D); sysd = c2d(sys, Ts);
    Ad = sysd.A; Bd = sysd.B; Cd = sysd.C; Dd = sysd.D; 
 
    y0 = 0.5; r0 = 0; x = [];
    x0 = [r0; y0];
    x(:,1) = x0;
    y(1) = Cd * x(:,1);

    D1 = zeros(Nu,Nu);

    temp = 1;
    for i=1:Nu
     F1(i,:) = Cd * Ad^(i-1);
    end
    temp = 1;
    for i=1:Nu
     D_(Nu-i+1) = Cd * temp * Bd;
     temp = temp * Ad;
    end
    for i=1:Nu
     D1(Nu-i+1,:) =  [D_(i+1:end) zeros(1,i)];
    end

    downsampling_factor = Fsu/Fsy;

    A = downsample(F1,downsampling_factor);
    B = downsample(D1,downsampling_factor);
end


%% stagnation handling
% function [tau_est,ub,lb]= handle_stagnation_in_boundary(y1,y2, sigma1, sigma2, u, tau, ub, lb, tol, tu, Fsy, Fsu)
%     tau_ = tau;
%     tau_est = tau;
%     J0 = Inf;
%     maxiter = 32;
%     range_of_search = 0.075;
%     ty = 0:(1/Fsy):(length(y1)-1)*(1/Fsy);
%     tu = 0:(1/Fsu):(length(u)-1)*(1/Fsu);
%     if(abs(tau(1)-ub(1))<tol || abs(tau(1)-lb(1))<tol)
%         for i=-maxiter:maxiter
%             b = tau(1);
%             tau_(1) = (b+((i)/maxiter)*range_of_search);
%             tau_=tau_(:);
%             J = cost_function_interior_point([tau_(1:2,:);tau_(1:2,:);tau_(3,:)],y1, y2, sigma1, sigma2, ty, u, tu);
%             if(J < J0 )
%                 tau_est = tau_;
%                 J0 = J;
%             end
%         end
%     end
%     if(abs(tau(2)-ub(2))<tol*3 || abs(tau(2)-lb(2))<tol*3)
%         for i=-maxiter:maxiter
%             b = tau(2);
%             tau_(2) = (b+((i)/maxiter)*range_of_search);
%             tau_=tau_(:);
%             J = cost_function_interior_point([tau_(1:2,:);tau_(1:2,:);tau_(3,:)],y1, y2, sigma1, sigma2, ty, u, tu);
%             if(J < J0 )
%                 tau_est = tau_;
%                 J0 = J;
%             end
%         end
%     end
%     ub(1:2) = [max(tau(1),ub(1)); max(tau(2),ub(2))];
%     lb(1:2) = [min(tau(1),lb(1)); min(tau(2),lb(2))];
% end
