%% comparative matrix
% Idx=(CR(1,:))'; %sIB�У�ÿ��ȡ���ǵ���������ÿһ�У������ÿ�ε����Ľ��
Idx=CR;%Idx�����վ�������ά��Ϊ1*N
D=gnd;%gnd����ʵ�����ݣ�ά��ΪN*1

%% AC evaluation      
 %CL = bestMap(D,Idx');
   AC = averageAC(Idx,D) ;
   
 %% MI evaluation
 num=size(CR,1);
 MIhat=zeros(1,num);
 for i=1:num
    MIhat(i)= MutualInfo(D,CR(i,:)); % pay attention to the CR(i,:);
 end
 MIhat=sum(MIhat,2)/num;
    
        
 %% ���ɻ������� 
 %conM = confusionmatrix(indx,L,clusters,categories)
 %clusters �Ǿ�������������category��ԭ���������
 %index�Ǿ�������L��data.L1���ǩ
 conM = confusionmatrix(Idx',D,2,2);
           
 %% ���㾫ȷ���������
 %Pre = sum(max(conM))/length(CL);
        