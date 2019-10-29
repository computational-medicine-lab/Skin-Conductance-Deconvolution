%codeForExperiment
% If you are using this code please cite the following paper

% Wickramasuriya DS, Amin MR and Faghih RT (2019) 
% Skin Conductance as a Viable Alternative for Closing the Deep Brain Stimulation Loop in Neuropsychiatric Disorders. 
% Front. Neurosci. 13:780. doi: 10.3389/fnins.2019.00780

close all;
clear all;
load('UT_Dallas_data\s.mat') % load UT dallas



sub_fron_neuro = [1, 5, 8, 9, 12, 16]; %Frontiers in Neuroscience subjects

%     Time        index                         Experimental Condition 
%----------------------------------------------------------------------
%     0:00.125        1     "    0    0    0	Relax
%     5:00.125     2401     "    0    0    0	PhysicalStress
%    10:26.250     5010     "    0    0    0	Relax
%    15:26.375     7411     "    0    0    0	EmotionalStress
%    16:06.500     7732     "    0    0    0	CognitiveStress
%    22:01.625    10573     "    0    0    0	Relax
%    27:01.750    12974     "    0    0    0	EmotionalStress
%    32:59.875    15839     "    0    0    0	Relax
downsampling_factor = 4;

ub = [1.4 6]';
lb = [0.1 1.5]';
dirname = 'deconvolution_results';
mkdir(dirname);

for subject = sub_fron_neuro(1)
    
    %% preprocessing
    % lowpass filter with 64 order 0.5 Hz filter
    Fs = 8;                    % sample rate in Hz
    N = length(s(subject).x);  % number of signal samples
    y = downsample(s(subject).x,downsampling_factor); % we downsample the signal to have a 2 Hz signal

    Fs = Fs/downsampling_factor;% new sampling frequency
    rng default;
    
    % Design a 64th order lowpass FIR filter with cutoff frequency of 75 Hz.
    Fnorm = 0.5/(Fs/2);           % Normalized frequency
    df = designfilt('lowpassfir','FilterOrder',64,'CutoffFrequency',Fnorm);
    D = mean(grpdelay(df)); % mean group delay
    temp = filter(df,[y'; zeros(D,1)]); % Append D zeros to the input data in order to adjust for group delay later
    temp = temp(D+1:end)';                  % Shift data to compensate for delay
    
    %% tonic component removal
    % [xn, mu, sigma] = zscore(temp);
    [r, p, t, l, d, e, obj] = cvxEDA(temp, 1/Fs, 2, 0.7, 10, 40e-4, 1e-3, 'quadprog');

    %% take 3*60 seconds of cognitive stress 
    ll = round(s(subject).y(5)/downsampling_factor);
    rr = round(s(subject).y(5)/downsampling_factor+3*60*Fs);
    Ts = 1/Fs;
    ty = 0:Ts:(length(ll:rr)-1)*Ts;
    CogStressMath.actual_signal = temp(ll:rr)';
    CogStressMath.tonic_part = t(ll:rr);
    CogStressMath.phasic_part = r(ll:rr);
    CogStressMath.phasic_driver = p(ll:rr);

    tic
    phasic_part = CogStressMath.phasic_part;
    Fsu = 4;
    Fsy = Fs;
    minimum_peak_distance = 1;  
    parallal_operations = 16;
    parfor i=1:parallal_operations
    [results(i).tau_j, results(i).uj, results(i).y, results(i).lambda, results(i).convergenceFlag] = coordinate_descent1(phasic_part, ub, lb, Fsu, Fsy, minimum_peak_distance);
    end
    toc
    cost_prev1 =Inf;
    cost_prev2 =Inf;
    y1 = phasic_part;
    
    % comparing all the initializations to take the best one
    for i=1:parallal_operations
        tau_j_ = results(i).tau_j;
        uj_ = results(i).uj;
        lambda_ = results(i).lambda;
        Nu = length(uj_);
        [A1, B1] = create_A_B_matrix_ss_multires(tau_j_(1:2), Nu, Fsu, Fsy);
        y = y1;
        y_ = A1*[0;y(1)] + B1 * uj_;
        
        % we consider two cost functions for comparison

        cost1 = 0.5 * norm(y-y_,2).^2;
        cost2 = 0.5 * norm(y-y_,2).^2 + lambda_ * norm(uj_, 1);

        if(cost1<cost_prev1 && results(i).convergenceFlag == 1 && round(results(i).tau_j(1)*1e4)/1e4 ~= lb(1) && round(results(i).tau_j(1)*1e4)/1e4 ~= ub(1))
            results1 = results(i); 
            cost_prev1 = cost1;
        end
        if(cost2<cost_prev2 && results(i).convergenceFlag == 1 && round(results(i).tau_j(1)*1e4)/1e4 ~= lb(1) && round(results(i).tau_j(1)*1e4)/1e4 ~= ub(1))
            results2 = results(i);
            cost_prev2 = cost2;
        end
    end
    
    % "result1" is taken as the solution
    
    ll = round(s(subject).y(4)/downsampling_factor);
    rr = floor(s(subject).y(end)/downsampling_factor+5*60*Fs)-2;
    y_w = r(ll:end);
    
    Nu = length(y_w)*Fsu/Fsy;
    
    % solving inverse problem only for neural stimuli (u(t)) for entire signal
    uj_1 = ones(Nu,1);
    [A1, B1] = create_A_B_matrix_ss_multires(results1.tau_j(1:2), Nu, Fsu, Fsy);
    y_w_ = y_w-A1*[0;y_w(1)];
    [uj_whole] = focuss_modified2(uj_1, y_w_, B1, -1, true, 15, 0.5, 10e-4, 1e-1 , round(minimum_peak_distance * Fsu));
    [uj_whole, J, Reg] = focussreg3_modified(uj_whole, y_w_, B1, 1e-2);

    save(strcat([dirname,'\result_s'],num2str(subject)));
end