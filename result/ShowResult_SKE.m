clc;clear;
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};
  
  %fprintf("  Name   |    RTR    |  HSODM-1   |  HSODM-10  |   HSODM-D2 \n")
  fprintf("  Name   |    RTR    |  HSODM  \n")

 for N = [1000,2000,5000]
    for j = [20,50]
    %name  = split(set{i},'.');
    name = "n"+num2str(N)+"r"+num2str(j);
    T     = load("manopt\SKE\"+name+"-result.mat");
    Out1  = load("hsodm\SKE\"+name+"-result.mat");
    fprintf("%6s      %3d       %3d\n",name{1},T.info(end).iter,Out1.Out.iter);
    end
 end