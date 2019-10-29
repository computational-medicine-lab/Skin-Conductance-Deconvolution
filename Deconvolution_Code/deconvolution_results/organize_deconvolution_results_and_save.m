
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Organize_Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;
fs = 2;
fsu = 4;
fss = 2;
bw = fsu;
subjects = [1, 5, 8, 9, 12, 16];
end_ = 5*60*fsu;
subject = [];

idx_count = 0;

for i = subjects
    idx_count = idx_count+1;
    data = load(['D:\rafi\SC_deconvolution_with_coordinate_descent\PlosOne\UT Dallas_with_deconvolution\results_12_04_2018_5\result_s',num2str(i),'.mat']);

    
    deconv_result_on_3min_segment.u = data.results1.uj;
    deconv_result_on_3min_segment.tau_r = data.results1.tau_j(1);
    deconv_result_on_3min_segment.tau_d = data.results1.tau_j(2);
    deconv_result_on_3min_segment.y_segment = data.results1.y;

    [A,B] = create_A_B_matrix_ss_multires(data.results1.tau_j(1:2), length(data.results1.uj), fsu, fs);
    deconv_result_on_3min_segment.y_reconstructed_segment = A*[0;data.results1.y(1)]+B*data.results1.uj;

    sparse_recovery_on_15min_40sec_segment.u = data.uj_whole(1:end-end_); sparse_recovery_on_15min_40sec_segment.u = sparse_recovery_on_15min_40sec_segment.u(:);
    
    sparse_recovery_on_15min_40sec_segment.y_segment = data.y_w;

    subject(idx_count).deconv_result_on_3min_segment = deconv_result_on_3min_segment;
    subject(idx_count).sparse_recovery_on_15min_40sec_segment = sparse_recovery_on_15min_40sec_segment;
    subject(idx_count).starting_time_of_tasks = data.s(idx_count).y;
    subject(idx_count).EDA_raw = data.s(idx_count).x;
    subject(idx_count).tonic_from_cvx_EDA = data.t;
    subject(idx_count).phasic_from_cvx_EDA = data.r;
    subject(idx_count).ID = i;
    
    save('Deconvolution_on_Experimental_data.mat','subject');
    
end