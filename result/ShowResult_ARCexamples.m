clc;clear;close all;
fprintf("  Name   |   HSODM   |   RTR  |   BFGS  |   BB  |   CG  \n")
Out_manopt =load("Result_manopt.mat");

datapath = "hsodm\";

width  = 10;     % Width in inches
height = 6;    % Height in inches
alw    = 0.75;    % AxesLineWidth
fsz    = 12;      % Fontsize
lw     = 2;      % LineWidth
msz    = 8;       % MarkerSize

Out_hsodm_DIS      = load(datapath+"Result_HSODM_DIS");
Out_hsodm_lrmc     = load(datapath+"Result_HSODM_lrmc");
Out_hsodm_rotation = load(datapath+"Result_HSODM_rotation");
Out_hsodm_shapefit = load(datapath+"Result_HSODM_shapefit");
Out_hsodm_SVD      = load(datapath+"Result_HSODM_SVD");
Out_hsodm_Maxcut   = load(datapath+"Result_HSODM_Maxcut");

Out_hsodm = {Out_hsodm_DIS,Out_hsodm_SVD,Out_hsodm_lrmc,Out_hsodm_Maxcut,Out_hsodm_rotation,...
            Out_hsodm_shapefit};


numproblems = 6;
numsolvers  = 6;
count = 1; 
for i = 1:numproblems
    subplot(3,2,count);
    for j = 1:2%numsolvers
        iter = length(Out_manopt.Outs{i,j}.gradnorm);
        semilogy(1:iter,Out_manopt.Outs{i,j}.gradnorm,'-o','LineWidth',1);
        %ylim([10^(-11),Out_manopt.Outs{i,j}.gradnorm(1)]);
        yticks([10^(-5),1]);
        %yticklabels({'10^(-10)','10^(-5)','10'})
        xlim([1,iter]);
        set(gca,'TickLabelInterpreter','latex' ,'FontSize', fsz, 'LineWidth', alw); %<- Set properties
        hold on;
    end
    semilogy(1:Out_hsodm{i}.Out.iter,Out_hsodm{i}.Out.grad,'-o','LineWidth',1);
    count = count +1;
    xlabel('iteration','Interpreter','latex');
    ylabel('$\|\nabla f(x_k)\|$','Interpreter','latex');
    switch i
        case 1
            title('Dominant invariant subspace','interpreter','latex','FontSize', fsz);
        case 2
            title('Truncated SVD','interpreter','latex','FontSize', fsz);
        case 3
            title('Low-rank matrix completion','interpreter','latex','FontSize', fsz);
        case 4
            title('Max-Cut','interpreter','latex','FontSize', fsz);
        case 5
            title('Synchronization of rotations','interpreter','latex','FontSize', fsz);
        case 6
            title('ShapeFit','interpreter','latex','FontSize', fsz);
    end
end
hold on




set(gcf, 'Position', [300 100  width*100, height*100]); %<- Set size
%legend('RTR','ARC','HSODM');

lg = legend( 'RTR',...
    'ARC', ...
    'RHSODM','interpreter','latex','NumColumns',4,'Position',[0.115,-0.01,0.8,0.1],'FontSize', fsz);

legend('boxoff');




%print("Numerical-result",'-depsc','-tiff');

hold off
count = 1;
figure();
for i = 1:numproblems
    subplot(3,2,count);
    for j = 3:numsolvers
        iter = length(Out_manopt.Outs{i,j}.gradnorm);
        semilogy(1:iter,Out_manopt.Outs{i,j}.gradnorm,'LineWidth',1);
        hold on;
    end
    %semilogy(1:Out_hsodm{i}.Out.iter,Out_hsodm{i}.Ou.grad);
    count = count +1;
    legend('CG','BB','GD','BFGS');
    xlabel('iteration','Interpreter','latex');
   % ylabel('\|\|','Interpreter','latex');
end


hold off
count = 1;
figure();
for i = 1:numproblems
    subplot(3,2,count);
    collectiters = [];
    for j = 1:numsolvers
        iter = length(Out_manopt.Outs{i,j}.gradnorm);
        collectiters = [collectiters,iter];
        semilogy(1:iter,Out_manopt.Outs{i,j}.gradnorm,'LineWidth',2);
        hold on;
    end
    semilogy(1:Out_hsodm{i}.Out.iter,Out_hsodm{i}.Out.grad,'LineWidth',2);
    %semilogy(1:Out_hsodm{i}.Out.iter,Out_hsodm{i}.Out.grad);
    count = count +1;
    xlim([1,max(collectiters)]);
    legend('RTR','ARC','CG','BB','GD','BFGS','HSODM');
    xlabel('iteration','Interpreter','latex');
    ylabel('gradnorm','Interpreter','latex');
end



for N =5000%[1000,2000,5000]
    for j =50%[20,50]
    name = "n"+num2str(N)+"r"+num2str(j);
    load("..\data\KSE\"+name{1}+".mat");
    
    Out_RTR     = load("manopt\RTR\SKE\"+name{1}+"-result.mat");
    Out_BFGS    = load("manopt\BFGS\SKE\"+name{1}+"-result.mat");
    Out_BB      = load("manopt\BB\SKE\"+name{1}+"-result.mat");
    Out_CG      = load("manopt\CG\SKE\"+name{1}+"-result.mat");
    

    %Out_GD     = load("manopt\GD\MaxCut\"+name{1}+"-result.mat");
    Out_hsodm   = load("hsodm\SKE\"+name{1}+"-result.mat");
    
    fprintf("%6s      %3d       %3d      %3d      %3d      %3d\n", ...
        name{1},Out_hsodm.Out.iter,Out_RTR.info(end).iter, ...
        Out_BFGS.info(end).iter,Out_BB.info(end).iter,Out_CG.info(end).iter);
    
    
    %subplot(1,2,1);

    NormRTR = zeros(length(Out_RTR.info),1);
    for j = 1 : length(Out_RTR.info)
        NormRTR(j) = Out_RTR.info(j).gradnorm;
    end
    semilogy(1:length(NormRTR),abs(NormRTR));

    hold on;
    NormBFGS = zeros(length(Out_BFGS.info),1);
    for j = 1 : length(Out_BFGS.info)
        NormBFGS(j) = Out_BFGS.info(j).gradnorm;
    end
    semilogy(1:length(NormBFGS),abs(NormBFGS));

    NormBB = zeros(length(Out_BB.info),1);
    for j = 1 : length(Out_BB.info)
        NormBB(j) = Out_BB.info(j).gradnorm;
    end
    semilogy(1:length(NormBB),abs(NormBB));

    NormCG = zeros(length(Out_CG.info),1);
    for j = 1 : length(Out_CG.info)
        NormCG(j) = Out_CG.info(j).gradnorm;
    end
    semilogy(1:length(NormCG),abs(NormCG));

    semilogy(1:length(Out_hsodm.Out.obj),Out_hsodm.Out.grad);
    
    xlabel('iteration','Interpreter','latex');
    ylabel('gradnorm','Interpreter','latex');
    legend('RTR','BFGS','BB','CG','HSODM');
%     
%     hold on;
% 
%     CostGD = zeros(length(Out_GD.info),1);
%     for j = 1 : length(Out_GD.info)
%         CostGD(j) = Out_GD.info(j).cost;
%     end
%     %semilogy(1:length(CostGD),abs(CostGD-Truecost));
% 
%     semilogy(1:length(Out_hsodm.Out.obj),Out_hsodm.Out.obj-Truecost);
%     %legend('RTR','GD','HSODM');
%     legend('RTR','HSODM');
%     title('Cost value gap');
%     subplot(1,2,2);
%     GradRTR = zeros(length(Out_RTR.info),1);
%     for j = 1 : length(Out_RTR.info)
%         GradRTR(j) = Out_RTR.info(j).gradnorm;
%     end
%     semilogy(1:length(GradRTR),abs(GradRTR));
% 
%     hold on;
% 
%     GradGD = zeros(length(Out_GD.info),1);
%     for j = 1 : length(Out_GD.info)
%         GradGD(j) = Out_GD.info(j).gradnorm;
%     end
%     %semilogy(1:length(GradGD),abs(GradGD));
% 
%     semilogy(1:length(Out_hsodm.Out.grad),Out_hsodm.Out.grad);
%     %legend('RTR','GD','HSODM');
%     legend('RTR','HSODM');
%     title('gradnorm');
    end
 end