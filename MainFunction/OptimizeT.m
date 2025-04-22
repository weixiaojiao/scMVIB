function  [T] = OptimizeT (T,cellX,cellInp,prm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

LoopCounter = 1;
NoChange =0;
done=0;

%% intial weights
Weights=ones(prm.m,1)*(1/(prm.m));

while ~done
    change=0;
    %% updating T
    for x=1:size(cellX{1},1)                % documents in every dataset should be same
        [tmp,t,inds_t]=ProcessTmp(cellInp,T,prm,x);
        if length(inds_t)>1
            [tmp] =Reduce_x (tmp,t,inds_t,x,cellInp,prm);  %change every Pt and P(y,t)
            [cellCosts] =MergeCosts (tmp,prm);
            [~,new_t] =MinSumCosts(cellCosts,prm,Weights);
            if new_t==t
                NoChange= NoChange+1;
            else
                NoChange= 0;
                change= change+1;
                [T] =UpdateAssignment (T,x,new_t,tmp,prm);
            end
        end  %end if length
    end %end for x=1:
    [T, done] =CheckConvergence (T,prm,LoopCounter,change);
    LoopCounter= LoopCounter+1;
    

    %% New~
    %% Update weights W:
    % computing objective
    ITY=zeros(prm.m,1);
    ITX=zeros(prm.m,1);
    HT=zeros(prm.m,1);
   % 初始化权重数组
   w = zeros(prm.m, 1);


for i=1:prm.m
        dd=T{i}.Pt ;
        PPt = repmat(dd,size(cellX{i},2),1);   % PPT is used to calculate P(t,y) in the next step
        Pty = T{i}.Py_t.*PPt;
	Ptx = T{i}.Pt_x.*PPt;
        [aa,~,cc] = MTC_local_MI (Ptx,dd);
        [bb,~,~] = MTC_local_MI (Pty,dd);
        ITX(i)=aa;
        ITY(i)=bb;
        HT(i)=cc;

	% 方法1：直接调用（假设 Ptx 列和已等于 P(x)）
    	%[ITX(i), ~, ~] = MTC_local_MI(Ptx, dd);
    
	% 方法2：显式对齐 Ptx 列和与 P(x)（更安全）
	% Ptx = Ptx ./ sum(Ptx, 1) .* cellInp{i}.Px;
	% [ITX(i), ~, ~] = MTC_local_MI(Ptx, dd);

    end


    % 遍历每一个数据点，计算权重
    for i = 1:prm.m
    	 % 根据公式计算权重
    	 w(i) = 1 / sqrt(abs(ITX(i) - prm.inv_beta * ITY(i)));
    end
end  %end while

% this also needed to be changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for i=1:prm.m
    dd=T{i}.Pt ;
    PPt = repmat(dd,size(cellX{i},2),1);   % PPT is used to calculate P(t,y) in the next step
    Pty = T{i}.Py_t.*PPt;
    [aa,bb,cc] = MTC_local_MI (Pty,dd);
    T{i}.Ity=aa;
    T{i}.Hy=bb;
    T{i}.Ht=cc;
    T{i}.L = T{i}.Ity;
%     T{i}.L = T{i}.Ity - prm.inv_beta*T{i}.Ht;
    % There seems to be a saying that I (T;X) = Ht
end

end
