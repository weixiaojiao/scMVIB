function conM = confusionmatrix(indx,L,clusters,categories)
%clusters �Ǿ�������������category��ԭ���������
%index�Ǿ�������L��data.L1���ǩ

conM = zeros(categories,clusters);
u = unique(indx);  
%ȡ����a�Ĳ��ظ� Ԫ�� ���ɵ�����
for i = 1:length(indx)  
     %length(x)��ʾ����x����������ĸ�������������ά��
     %������������������֮����ֵ�ϴ����һ��ֵ�����൱��ִ����max(size(A))
    r = find(L(i,:) ~= 0);
    %MATLAB��"��"��ʹ�÷ǳ����򵥵�˵���Ǵ�������е�ȫ�л�ȫ��Ԫ��
    %ð�ű�ʾ���е�index�����Y(i,j)��ʾY�ĵ�i�е�j�У�Y(:,i)�ͱ�ʾ������i����
    c = find(u == indx(i));
    %index�е�i ��Ԫ��
    %����find�����������������Ԫ��λ��
    conM(r,c) = conM(r,c) + 1;
end