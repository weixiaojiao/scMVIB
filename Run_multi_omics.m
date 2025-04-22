clear all
path1='MainFunction';
path2='BasicFunction';
path3='clustering performance';
addpath(path1,path2,path3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auto-weighted multi-view IB method.
dataset_str = 'SLN206D1';
m = 3;               % num of the input view
Restarts=1;

Path = 'D:/Tesedemo/from scMDC/CITEseq_sln206/Batch01/';
rna_path = strcat(Path,'rna_cite_sln206_Batch1.tsv');
label_path = strcat(Path,'label_cite_sln206_Batch1.tsv');
view2_path = strcat(Path,'adt_cite_sln206_Batch1.tsv');
IE_path = append('D:/Experiment/Code/Python Code/scRNAseq code package/MOINER-master/results_preprocessing/',dataset_str,'/4.Multi-Omics_Data(Information Enhanced).csv');

Select_beta = 50

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellX = cell(m,1);     % input data
UniformPrior=1;
CentroidFigs=0;
Select_theta=5

% 构造input data
expr_matrix = importfile_tsv(rna_path);
sum_values = sum(expr_matrix,2);
index = find(sum_values ~= 0);
expr_matrix = expr_matrix(index,:);
expr_matrix = sparse(expr_matrix');
cellX{1} = expr_matrix;

% ATAC_data = importdata(view2_path);  % delete this 

view2_data = importfile_tsv(view2_path);
sum_values = sum(view2_data,2);
index = find(sum_values ~= 0);
view2_data = view2_data(index,:);
view2_matrix = view2_data;
view2_matrix = sparse(view2_matrix');
cellX{2} = view2_matrix;


% IE matrix:
IE_matrix = importfile_csv(IE_path);
IE_matrix = IE_matrix';  % IE matrix 原始：row：cell colnum：gene
sum_values = sum(IE_matrix,2);
index = find(sum_values ~= 0);
IE_matrix = IE_matrix(index,:);
IE_matrix = sparse(IE_matrix');
cellX{3} = IE_matrix;


% label files:
label = import_label(label_path);  % load label
k_num = length(unique(label));  % cluster num
label = grp2idx(label);  % label: char to double % 字符型标签转数值型标签
gnd = label;

% Run IB_Multi_version
for p=1:length(Select_theta)
    diary on;
    for q=1:length(Select_beta)
        tic;  % Running time
        [CR,cellPI,cellT,seed]=MainIB (cellX, m,k_num,Select_theta(p),Select_beta(q), Restarts, UniformPrior, CentroidFigs);
        toc;
        for j=1:Restarts   % record the result of each Restart, One Restart is equivalent to running one IB
            NMI(j)=func_nmi(gnd',CR(j,:));
            ARI(j)=func_ari(gnd',CR(j,:));
        end
        MeanARI=mean(ARI);
        StdARI=std(ARI);
        MeanNMI=mean(NMI);
        StdNMI=std(NMI);
        nmi = func_nmi(gnd',cellT{1}.Pt_x);
        ari = func_ari(gnd',cellT{1}.Pt_x);

        % the best result
        % % prepare data
        [m_ARI,p_ARI] = max(ARI);
        [m_NMI,p_NMI] = max(NMI);
        if p_ARI == p_NMI
            result = [m_ARI,m_NMI];
        else
            result = [m_ARI,NMI(p_ARI),ARI(p_NMI),m_NMI];
        end
        % result
        disp([' the view num is: ' num2str(m)]);
        disp([' the theta is: ' num2str(Select_theta(p))]);
        disp([' the Select_beta is: ' num2str(Select_beta(q))]);
        disp([' the Dataset and Restarts is:  ' dataset_str '  ' num2str(Restarts)]);
        disp([' the best result is: ' num2str(result)]);
        disp([' the final best ARI & NMI based on L is:  ' num2str(ari) '  ' num2str(nmi)]);
        disp([' ARI mean and std is:  ' num2str(MeanARI) '  ' num2str(StdARI)]);
        disp([' NMI mean and std is:  ' num2str(MeanNMI) '  ' num2str(StdNMI)]);
        disp([' Seed is:  ' num2str(seed)]);

        disp(['运行时间：' num2str(toc)]);
        fprintf('运行时间：%d小时 %d分钟 %d秒\n', floor(toc/3600),floor(mod(toc,3600)/60),mod(mode(toc,3600),60));

        save_path = append('Result/multi_omics/',dataset_str,'/',num2str(Select_theta(p)),'_',num2str(Select_beta(q)),'_',num2str(Restarts),'_', num2str(m),'_CR.csv');
        writematrix(CR(p_ARI,:)',save_path) % param :α β restart


    end
    disp('-----------分割线----------');
    diary off;
end

clear path1 path2 path3
clear data_path top2000_path KNN_path label_path
clear I index label sum_values
clear j q CentroidFigs UniformPrior
clear p_NMI p_ARI