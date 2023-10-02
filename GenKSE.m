clc; clear all;
rng('default');
for n = [1000,2000,5000,7000,10000]
    for r = [20,50]
        rng('default');
        %% generate data
        e = ones(n,1);
        L = spdiags([-e,2*e,-e],[-1,0,1],n,n);     
        save(["data\\KSE\\"+"n"+ num2str(n)+"r"+ num2str(r)]);
    end
end
