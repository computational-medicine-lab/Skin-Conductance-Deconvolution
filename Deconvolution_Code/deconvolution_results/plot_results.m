%plot results from experiment 2_8_18
% 4 Hz with support of 0.5 sec
% lb = [0.65 2];
% ub = [0.85 4];
clc;
clear all;
close all;
addpath('D:\rafi\DrosteEffect-BrewerMap-b6a6efc');
addpath('..')
count = 0;
sub = 0;
c = brewermap(9,'Accent');
mkdir('FigDeconv')
for sub = [1,5,8,9,12,16]
    
    figure('units','normalized','outerposition',[0 0 2/3*0.75 2/3*0.75]),
    sub
    count = count+1;
    data = load(strcat('result_s',num2str(sub)));
    y_original = downsample(data.s(sub).x,data.downsampling_factor);
    results1 = data.results1;Fsy = data.Fsy; Fsu = data.Fsu; y = data.results1.y;
    subplot(212),
    amp_correction = 1;
    u = results1.uj;Nu = length(u);
    sum(u>0)
    u(u==0)=NaN;
    ty = (0:1:(length(results1.y)-1))/Fsy;
    %subplot(3,2,count), 
%     grid on;
    yyaxis left,
    plot(ty(:),results1.y(:),'*','Linewidth',1,'color',c(1,:)),
    [A,B] = create_A_B_matrix_ss_multires(results1.tau_j(1:2), Nu, Fsu, Fsy);
    set(gca,'fontsize',20)
    y_est = A*[0;y(1)]+B*results1.uj;
    tu = (0:length(u)-1)/Fsu;
    hold on, plot(ty(:),y_est(:),'linewidth',2,'color',0.5*c(2,:)),
    ylabel('SC (\muS)','fontsize',20)
    %ylim([-0.2 2.1]);
    fig = gca;
    fig.YColor = 0.8*c(1,:);
    R_2 = 1-var(results1.y(:)-y_est(:))/var(y)
    hold on, yyaxis right,ylabel('Amplitude','fontsize',20)
    h = stem(tu(:),u/amp_correction,'fill','linewidth',2,'color',c(5,:));%set(h, 'Marker', 'none')
    h = get(gca,'Children');
    %plot(tu(:),u/amp_correction,'MarkerSize',30)

    %title(['Example of Deconvolution on Experimental Data'])
    xlabel('Time (seconds)','fontsize',20)
    %ylabel('Amplitude','fontsize',20)
    fig = gca;
    fig.YColor = c(5,:);
    box on;
    xlim([0 180]); ylim([0 max(u)*1.2]); 
    
    subplot(211),
    ll = round(data.s(sub).y(5)/data.downsampling_factor);
    rr = round((data.s(sub).y(5)/data.downsampling_factor+3*60*data.Fs));
    plot(ty(:),y_original(ll:rr),'Linewidth',5,'color',1.2*c(1,:)), hold on;
    plot(ty(:),data.CogStressMath.tonic_part+data.CogStressMath.phasic_part,'Linewidth',1.5,'color',.7*c(1,:)), hold on;plot(ty(:),data.CogStressMath.tonic_part,'Linewidth',1.5,'color',1.2*c(7,:))
    ddelta = (max(data.CogStressMath.tonic_part) - min(data.CogStressMath.tonic_part));
    ylim([min(data.CogStressMath.tonic_part)-ddelta*0.1 min(data.CogStressMath.tonic_part)+ddelta*2]);
    xlim([0 180]);
    xlabel('Time (seconds)','fontsize',20),ylabel('SC (\muS)','fontsize',20), fig = gca;
    fig.YColor = 0.8*c(1,:);
    fig.XAxis.FontSize = 20;
    fig.YAxis.FontSize = 20;
    saveas(gcf,['FigDeconv\DeconvExpP',num2str(count),'.png']);
%     saveas(gcf,['FigDeconv\DeconvExpP',num2str(count),'.eps']);
    title(['Participant ',num2str(count)],'fontsize',20);
    
end
