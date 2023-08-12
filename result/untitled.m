clc;clear;
set= {'mcp100.mat','mcp124-1.mat','mcp124-2.mat','mcp124-3.mat','mcp124-4.mat','mcp250-1.mat','mcp250-2.mat','mcp250-3.mat','mcp250-4.mat',...
      'mcp500-1.mat','mcp500-2.mat','mcp500-3.mat','mcp500-4.mat'};
  
  fprintf("  Name   |    RTR    |  HDODM-1   |  HDODM-10 \n")
  for i = 1:length(set)
    name  = split(set{i},'.');
    T     = load("manopt\"+name{1}+"-result.mat");
    %Out1  = load("hsodm\eta1\"+name{1}+"-result.mat");
    %Out10 = load("hsodm\eta10\"+name{1}+"-result.mat");
    %fprintf("%6s      %6d       %6d        %6d \n",name{1},length(T.info),Out1.Out.iter,Out10.Out.iter);
    info    = T.info(length(T.info));
    xcost   = T.xcost;
    options = T.options;
    save("manopt\temp\"+name{1}+"-result.mat",'info','xcost','options');
  end