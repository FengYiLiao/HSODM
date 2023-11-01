clc;clear;close all;
%set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
%      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};
  
  %fprintf("  Name   |    RTR    |  HSODM-1   |  HSODM-10  |   HSODM-D2 \n")
  fprintf("  Name   |    RTR    |  HSODM   |    GD \n")
 for N = 1000%[1000,2000,5000]
    for j = 20%[20,50]
    %name  = split(set{i},'.');
    name = "n"+num2str(N)+"r"+num2str(j);
    %name        = split(set{i},'.');
    load("..\data\KSE\"+name{1}+".mat");
    
    Out_RTR     = load("manopt\RTR\SKE\"+name{1}+"-result.mat");
    %Out_GD     = load("manopt\GD\MaxCut\"+name{1}+"-result.mat");
    Out_hsodm   = load("hsodm\SKE\"+name{1}+"-result.mat");
    
    fprintf("%6s      %3d      %3d\n", ...
        name{1},Out_RTR.info(end).iter,Out_hsodm.Out.iter);
    
    
    subplot(1,2,1);

    CostRTR = zeros(length(Out_RTR.info),1);
    for j = 1 : length(Out_RTR.info)
        CostRTR(j) = Out_RTR.info(j).cost;
    end
    semilogy(1:length(CostRTR),abs(CostRTR-Truecost));
    
    hold on;

    CostGD = zeros(length(Out_GD.info),1);
    for j = 1 : length(Out_GD.info)
        CostGD(j) = Out_GD.info(j).cost;
    end
    %semilogy(1:length(CostGD),abs(CostGD-Truecost));

    semilogy(1:length(Out_hsodm.Out.obj),Out_hsodm.Out.obj-Truecost);
    %legend('RTR','GD','HSODM');
    legend('RTR','HSODM');
    title('Cost value gap');
    subplot(1,2,2);
    GradRTR = zeros(length(Out_RTR.info),1);
    for j = 1 : length(Out_RTR.info)
        GradRTR(j) = Out_RTR.info(j).gradnorm;
    end
    semilogy(1:length(GradRTR),abs(GradRTR));

    hold on;

    GradGD = zeros(length(Out_GD.info),1);
    for j = 1 : length(Out_GD.info)
        GradGD(j) = Out_GD.info(j).gradnorm;
    end
    %semilogy(1:length(GradGD),abs(GradGD));

    semilogy(1:length(Out_hsodm.Out.grad),Out_hsodm.Out.grad);
    %legend('RTR','GD','HSODM');
    legend('RTR','HSODM');
    title('gradnorm');
    end
 end