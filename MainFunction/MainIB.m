function [CR,PI,T,seed] = MainIB(cellX, m,Tsize,theta, beta, Restarts, UniformPrior, CentroidFigs)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

format compact

%init function
[prm]=Initparameters(m,Tsize,theta, beta, Restarts, UniformPrior, CentroidFigs);
[cellInp]=ProcessInput(cellX,prm);


%main loop
CR=zeros(Restarts,size(cellX{1},1)); %clustering result
PI=zeros(1,Restarts);      %final optimal L
bestL=-inf;

for RS=1:prm.Restarts
%     fprintf ('Restart number %d....\n',RS);
    [tmpT]= RandomPartition(cellX,cellInp,prm);
    [tmpT] = OptimizeT (tmpT,cellX,cellInp,prm);
    
    CR(RS,:)=tmpT{1}.Pt_x;
    PI(RS)=tmpT{1}.L(end);
    if tmpT{1}.L(end) > bestL
        T{1}= tmpT{1};
        bestL= T{1}.L(end);             % use the final as the best,cause it is convergent
    end
end

% produce centroid figures (optional)

%  Last updates
% for i=1:prm.m
% cellT{i}.Ity_div_Ixy = cellT{i}.Ity(end)/cellInp{i}.Ixy;
% cellT{i}.Ht_div_Hx = cellT{i}.Ht(end)/cellInp{i}.Hx;
% end

seed = prm.RunSeed;

end

