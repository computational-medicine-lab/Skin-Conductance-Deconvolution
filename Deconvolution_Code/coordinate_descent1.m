function [tau_j, uj, y, lambda, convergenceFlag] = coordinate_descent1(y, ub, lb, Fsu, Fsy, minimum_peak_distance)

%assumption
%theta1 = theta2, alpha ~= 1;

%% approximate number of peaks
[~,locs]=findpeaks(y, 'MinPeakDistance',minimum_peak_distance * Fsy);
ubn=length(locs)+20;
%% random initialization of the system parameters
tau_j_1(1) = lb(1) + (ub(1)-lb(1))*rand(1,1);   tau_j_1(2) = lb(2) + (ub(2)-lb(2))*rand(1,1);

Nu = length(y) * Fsu/Fsy;
y0 = y(1);

Tsy = 1/Fsy; %min
Tsu = 1/Fsu; %min
ty = 0:Tsy:(length(y)-1)*Tsy;
tu = 0:Tsu:(Nu-1)*Tsu;

%% initialization 
for j=1:30
    
    [A, B] = create_A_B_matrix_ss_multires(tau_j_1, Nu, Fsu, Fsy);  % create dictionary matrices
    
    y_ = y-A*[0;y0]; % subtract the initial condition
   
    uj_1 = ones(Nu,1);  % initialize vector u with all ones
    [ uj ] = focuss_modified2(uj_1, y_, B, ubn, true, 15, 0.5, 1e-4, 1e-5 , round(minimum_peak_distance * Fsu)); % solve the inverse problem with sparse recovery to find the neural stimuli
    
    tau_j = SystemID_interior_point1(tau_j_1, y, ty, uj, tu, ub, lb); % estimate the system parameters
    tau_j_1 = tau_j;
end

%% coordinated descent
count = 0;
maxiter = 500; % set maximum iteration.
while(1)

    [A, B] = create_A_B_matrix_ss_multires(tau_j_1, Nu, Fsu, Fsy); % create dictionary matrices
    y_ = y-A*[0;y0]; % subtract the initial condition
    uj_1 = uj; % initialize vector with previous solution
    [uj, J, Reg] = focussreg3_modified(uj_1, y_, B, 0.01); % solve the inverse problem with generalized cross-validation based sparse recovery to find the neural stimuli
    tau_j = SystemID_interior_point1(tau_j_1, y, ty, uj, tu, ub, lb);
    tau_j_1 = tau_j;
    
    % display the estimated system arameters
%     fprintf('count = %d---->   ',count);
%     fprintf('tau1 = %d,   tau2 = %d\n', tau_j(1),tau_j(2));
    
    % check for convergence
    if(round((uj)*5e1)/5e1 == round((uj_1)*5e1)/5e1)
                disp('convergence achieved for u');
                 if(round(tau_j(1)*1e2)/1e2 == round(tau_j_1(1)*1e2)/1e2 && round(tau_j(2)*1e2)/1e2 == round(tau_j_1(2)*1e2)/1e2 && round(tau_j(2)*1e2)/1e2 == round(tau_j_1(2)*1e2)/1e2)
                     disp('convergence achived for tau');
                     convergenceFlag = 1;
                     break;
                 end
    end
    count = count+1;
    if(count>maxiter)
        convergenceFlag = 0;
        break;
    end
end

lambda = Reg;
end

%% system identification
function tau_est = SystemID_interior_point1(tau_prev, y, ty, u, tu, ub, lb)
    lb = lb(:); ub = ub(:);
    fun = @(x)cost_function_interior_point(x,y, ty, u, tu);
    %x0 = (ub+lb)/2;  
    x0 = tau_prev(:);
    A = []; b = [];  Aeq = [];  beq = [];
    nonlcon = [];
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off'); % run interior-point algorithm
    tau_est = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    tau_est = tau_est(:)';
end

%% cost function
function J = cost_function_interior_point(tau,y,ty, u, tu)
    Fsu = 1/(tu(2)-tu(1)); % retrieve neural stimuli frequency
    Fsy = 1/(ty(2)-ty(1)); % retrieve signal frequency
    tau = tau(:)';
    [A,B] = create_A_B_matrix_ss_multires(tau, length(tu), Fsu, Fsy); 
    y0 = [0;y(1)];
    J = 0.5*norm(y-A*y0-B*u,2)^2;
end

function [A, B] = create_A_B_matrix_ss_multires(tau, Nu, Fsu, Fsy)

% this function is created based on the state-space approach
% the detailed description of the theory can be found in the following
% papers (see supplemnetary materials)

% Faghih, R.T., Dahleh, M.A., Adler, G.K., Klerman, E.B. and Brown, E.N., 2015. 
% Quantifying pituitary-adrenal dynamics and deconvolution of concurrent cortisol and adrenocorticotropic hormone data by compressed sensing. 
% IEEE Transactions on Biomedical Engineering, 62(10), pp.2379-2388.

% Amin, M.R. and Faghih, R.T., 2019. 
% Sparse Deconvolution of Electrodermal Activity via Continuous-Time System Identification. 
% IEEE Transactions on Biomedical Engineering.


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