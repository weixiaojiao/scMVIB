function AC = averageAC(CR,gnd)  

%% �Զ��־�������ƽ��ֵ���Լ�С���,ͬʱҲ���÷�����⼸��ƽ��ֵ��Ĳ�����С��
%CR��ʾ��������gnd��ʾ��ʵ��Դ����;
T = size(CR,1);
AC = zeros(1,T);

for i = 1:T
    res = bestMap(gnd,CR(i,:));
    AC(i) = length(find(gnd == res))/length(gnd);
end

%% get the average
%  AC=sum(AC,2)/T;


    