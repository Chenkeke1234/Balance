clear
AvLoss3=zeros(1,5);%���㶪ʧ���ݰ�
AvLoss4=zeros(1,5);
AvLoss5=zeros(1,5);
AvLoss6=zeros(1,5);
AvLoss7=zeros(1,5);
VaLoss3=zeros(1,5);%���㶪ʧ���ݰ�
VaLoss4=zeros(1,5);
VaLoss5=zeros(1,5);
VaLoss6=zeros(1,5);
VaLoss7=zeros(1,5);
for o=1:1:5
PL=5.8;
AS3=zeros(1,100);%ͳ��S3��Ծ�ڵ�
Loss3=zeros(1,100);%���㶪ʧ���ݰ�
Loss4=zeros(1,100);
Loss5=zeros(1,100);
Loss6=zeros(1,100);
Loss7=zeros(1,100);
for q=1:1:100
    q
%1.��ʼ�����趨ģ�� 
%.�������ڵ��������(��λ M)
xm=100;
ym=100;
B=100;
rmax=1;
%(1)��۽��������
PI=[];
PI=0:pi/150:100*pi;
for i=1:1:rmax+1
    Sink(i).x=xm*cos(PI(i))*rand(1,1);
    Sink(i).y=ym*sin(PI(i))*rand(1,1);
end
sink.x=0;
sink.y=0;
%�����ڴ�������
n=100;
%��ͷ�Ż���������ѡ��ͷ�ĸ��ʣ�
P=0.05;
%����ģ�ͣ���λ ����
%��ʼ������ģ��
Eo=0.1;
%Eelec=Etx=Erx
ETX=50*0.000000001;%����������ÿbit
ERX=50*0.000000001;%����������ÿbit
%Transmit Amplifier types
Efs=10*0.000000000001;%���ɿռ���ģ�ÿbit
Emp=0.0013*0.000000000001;%�ྶ˥����ģ�ÿbit
%Data Aggregation Energy
EDA=5*0.000000001;%�ں��ܺģ�ÿbit
EPC=5*0.000000001;%ÿ�ֵĸ�֪�ռ���λ��Ϣ���ܺģ�ÿbit
%���ѭ������

%������� do���ڵ�ͨ�Ű뾶
do=sqrt(Efs/Emp)/2;
fo=ceil(do);
Et=0;
%���ݴ������
data3=0;
data4=0;
data5=0;
data6=0;
data7=0;
%2.���ߴ���������ģ�Ͳ���ģ��
%�������ߴ���������,�������ھ���Ͷ��100���ڵ�,������ͼ��
%��ʼ���ܼ�����
e=zeros(4,B);%4��ʾ������
a=zeros(4,B+1);
eav=zeros(4,B+1);
for i=1:1:n
    S3(i).xd=(rand(1,1)*2-1)*i;
    S3(i).yd=(rand(1,1)*2-1)*i;
    S4(i).xd=S3(i).xd;
    S4(i).yd=S3(i).yd;
    S5(i).xd=S3(i).xd;
    S5(i).yd=S3(i).yd;
    S6(i).xd=S3(i).xd;
    S6(i).yd=S3(i).yd;
    S7(i).xd=S3(i).xd;
    S7(i).yd=S3(i).yd;
    %����ÿ���ڵ����ڵ�ľ���
    S3(i).D=sqrt(S3(i).xd^2+S3(i).yd^2);
    S4(i).D=S3(i).D;
    S5(i).D=S3(i).D;
    S6(i).D=S3(i).D;
    S3(i).C=0;
    S6(i).C=0;
    S4(i).G=0;    
    S5(i).G=0;
    S6(i).G=0;
    S7(i).G=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if (S3(i).yd>=S3(i).xd && S3(i).yd>(-S3(i).xd))
         S3(i).C=1;
         S6(i).C=1;
    elseif (S3(i).yd<S3(i).xd && S3(i).yd>=(-S3(i).xd))
         S3(i).C=2;
         S6(i).C=2;
    elseif (S3(i).yd<=S3(i).xd && S3(i).yd<(-S3(i).xd))
         S3(i).C=3;
         S6(i).C=3;
    else
         S3(i).C=4;
         S6(i).C=4;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
     S3(i).f=0;
   if((abs(S3(i).xd)-abs(S3(i).yd))>=0)
     S3(i).f=ceil(abs(S3(i).xd));
   else
     S3(i).f=ceil(abs(S3(i).yd));
   end
    %�����ʼ����E0�Ǿ���һ�ַ��ͽڵ���Ϣ����ڵ��ʣ�������
    if(S3(i).D<=do)
        S3(i).E=Eo*(4-sqrt((S3(i).xd*S3(i).xd+S3(i).yd*S3(i).yd))/(xm*sqrt(2)));
        S3(i).EB=S3(i).E;
    else
        S3(i).E=Eo*(16-sqrt((S3(i).xd*S3(i).xd+S3(i).yd*S3(i).yd))/(xm*sqrt(2)));
        S3(i).EB=S3(i).E;
    end
    if(S3(i).D<=do) 
        S4(i).E=S3(i).E;
        S4(i).EB=S3(i).EB;
        S5(i).E=S3(i).E;
        S5(i).EB=S3(i).EB;
        S6(i).E=S3(i).E;
        S6(i).EB=S3(i).EB;
        S7(i).E=S3(i).E;
        S7(i).EB=S3(i).EB;
    else
        S4(i).E=sqrt(S3(i).E);
        S4(i).EB=sqrt(S3(i).EB);
        S5(i).E=sqrt(S3(i).E);
        S5(i).EB=sqrt(S3(i).EB);
        S6(i).E=sqrt(S3(i).E);
        S6(i).EB=sqrt(S3(i).EB);
        S7(i).E=sqrt(S3(i).E);
        S7(i).EB=sqrt(S3(i).EB);
    end
    E3(i)= S3(i).E;
    Et=Et+E3(i);%������
    %initially there are no cluster heads only nodes
    S3(i).type='N';%���c��͞���ͨ
    S4(i).type='N';
    S5(i).type='N';
    S6(i).type='N';
    S7(i).type='N';
    S3(i).M=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ھӽڵ��
L=10;%�ڵ��֪�뾶
R=50;%ͨ�Ű뾶
for i=1:1:n
    Contigupus=zeros(1,n);%�洢�ھӽڵ���
    Similar=zeros(1,n);
    A=zeros(1,n);
    for j=1:1:n
        d=sqrt((S3(i).xd-(S3(j).xd))^2 + (S3(i).yd-(S3(j).yd))^2);
        if(d<=2*L && d~=0 && S3(j).E>0)%�ڵ��������ڸ�֪�뾶��Ϊ�ھӽڵ�
            Contigupus(j)=j;   
        end
        if(d<L && d~=0 && S3(j).E>0)%��ѡ���ƽڵ�
            Similar(j)=j;   
        end
        if(S3(j).D>do && d<=do && d~=0 && S3(j).E>0 && S3(i).C==S3(j).C)
            A(j)=j;
        end
    end
    if (sum(find(Contigupus~=0))~=0)%�жϾ����Ƿ�Ϊ��
         S3(i).NB=find(Contigupus~=0);
         S7(i).NB=find(Contigupus~=0);
    else
         S3(i).NB='null';
         S7(i).NB='null';
    end
    if (sum(find(Similar~=0))~=0)%�жϾ����Ƿ�Ϊ��
         S3(i).S=find(Similar~=0);
    else
         S3(i).S='null';
    end
    if (sum(find(A~=0))~=0)%�жϾ����Ƿ�Ϊ��
         S3(i).a=length(find(A~=0))+1;
    else
         S3(i).a=1;
    end
end
d1=0.765*xm/2;
K=sqrt(0.5*n*do/pi)*xm/d1^2;
d2=xm/sqrt(2*pi*K);
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);

flag_first_dead3=0;
flag_teenth_dead3=0;
flag_all_dead3=0;
dead3=0;
first_dead3=0;
teenth_dead3=0;
all_dead3=0;
allive3=n;
%�����ڵ���
allive4=n;
flag_first_dead4=0;
flag_teenth_dead4=0;
flag_all_dead4=0;
dead4=0;
first_dead4=0;
teenth_dead4=0;
all_dead4=0;
%�����ڵ���
allive5=n;
flag_first_dead5=0;
flag_teenth_dead5=0;
flag_all_dead5=0;
dead5=0;
first_dead5=0;
teenth_dead5=0;
all_dead5=0;
%�����ڵ���
allive6=n;
flag_first_dead6=0;
flag_teenth_dead6=0;
flag_all_dead6=0;
dead6=0;
first_dead6=0;
teenth_dead6=0;
all_dead6=0;
%�����ڵ���
allive7=n;
flag_first_dead7=0;
flag_teenth_dead7=0;
flag_all_dead7=0;
dead7=0;
first_dead7=0;
teenth_dead7=0;
all_dead7=0;
 %%%%%%%%%%%%%%%%%%%%%
   %����С�㼶����
   for j=1:1:n
       e(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)+S3(j).E;%����С�㼶������
       a(S3(j).C,S3(j).f)=a(S3(j).C,S3(j).f)+1;%���������С�㼶�Ľڵ���
   end
   for j=1:1:n
       eav(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)./a(S3(j).C,S3(j).f);
       a(S3(j).C,B+1)= a(S3(j).C,B+1)+1;%���������Ľڵ����
   end
Ea=Et/n;%EaΪʣ����ƽ������
SCHCount=ceil(0.04*n);
N=zeros(SCHCount,n);%���澲̬��ͷɸѡ����
R=zeros(1,n);%���涯̬��ͷɸѡ����
JE=zeros(1,SCHCount);%�����ж������Ƿ���ڶ�̬��ͷ
countCHs3=0;%�����ͷ����ɸѡ��������ͷ��������
for i=1:1:n
 N(S3(i).C,i)=(eav(S3(i).C,S3(i).f)*(S3(i).E-Ea))*(S3(i).a/n)/(4*Eo-S3(i).E)*(1-(S3(i).D/(1.42*B)));
%  if Ea>0
%      p(i)=P*S3(i).a/n;
%      if(S3(i).E>0)    
%            %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
%         N(S3(i).C,i)=(p(i)/(1-p(i)*mod(0,round(1/p(i))))*(S3(i).E/(3*eav(S3(i).C,S3(i).f))));%p(i)<1,�����
%      end  
%   end 
end
[st,Sid]=max(N,[],2);%��ȡ��̬��ͷɸѡ�����е����ֵ����ID
for j=1:1:SCHCount
    if (st(j)~=0 && a(S3(Sid(j)).C,B+1)>10) %��������������������̬��ͷ
       S3(Sid(j)).type='C';
       S3(Sid(j)).Cid='sink';%û�ж�̬��ͷʱ����̬��ͷ���ڵ�sink����
       countCHs3=countCHs3+1;
       if( S3(Sid(j)).f>fo )
           for z=1:1:n 
               if(S3(z).D<S3(Sid(j)).D && S3(z).C==S3(Sid(j)).C)%ͬһ����͵���ڵ�������ھ�̬��ͷʱ���㶯̬�ڵ�Ӱ������
                   R(z)=(4*(S3(Sid(j)).f-S3(z).f+1)*S3(z).f/((S3(Sid(j)).f+1)*(S3(Sid(j)).f+1)))*(1-S3(z).f/S3(Sid(j)).f)*(eav(S3(z).C,S3(z).f)*(S3(j).E-Ea))/(3*Eo-S3(z).E);
               end
           end
       [dy,Did]=max(R,[],2); %��ȡ��̬��ͷɸѡ�����е����ֵ����ID  
       S3(Did).type='D';
       S3(Sid(j)).Cid=Did;%�ж�̬��ͷ������Ϊ�붯̬��ͷ����
       S3(Did).Cid='sink';%��̬��ͷ���ڵ�����
       JE(S3(Did).C)=Did;
       countCHs3=countCHs3+1;
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%�ڵ�����
LB=zeros(1,n);%������Ϣ������
for i=1:1:n
    if ( S3(i).type=='N' && S3(i).E>0 && S3(i).D<=do)
        S3(i).Cid='sink';
    elseif(S3(i).type=='N' && S3(i).E>0 && S3(i).D>do)
            if(JE(S3(i).C)~=0  && Sid(S3(i).C)~=0)%���ڶ�����̬��ͷ
                d1=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                d2=sqrt((S3(i).xd-S3(JE(S3(i).C)).xd)^2+(S3(i).yd-S3(JE(S3(i).C)).yd)^2);
                if(d1<d2)
                    S3(i).Cid=Sid(S3(i).C);%����ڵ������͵�����ڵ�ID
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                    S3(i).Cid=JE(S3(i).C);
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                end
            elseif(Sid(S3(i).C)~=0 && JE(S3(i).C)==0)%�����ڶ�̬��ͷ�����ھ�̬��ͷ
                d3=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                if(sum(S3(i).NB)==443 || d3<=do)%���û���ھӽڵ㣬��ֱ���뾲̬��ͷ����,sum('null')==443
                     S3(i).Cid=Sid(S3(i).C);
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                elseif(sum(S3(i).NB)~=443 && d3>do)
                    b=[];
                    for z=1:1:length(S3(i).NB)
                        b(z)=S3(S3(i).NB(z)).E;
                    end
                    for j=1:1:length(S3(i).NB)
                        if(S3(S3(i).NB(j)).type~='N')%���ھӽڵ���������ڵ�ʱ,��������ڵ�
                            S3(i).Cid=S3(i).NB(j);
                            LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            break;
                        else
                            if(S3(i).f>S3(S3(i).NB(j)).f && S3(S3(i).NB(j)).E==max(b))
                                 S3(i).Cid=S3(i).NB(j);%ѡ����㼶�������������������ھӽڵ�����
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                                 break;
                            else
                                 S3(i).Cid=Sid(S3(i).C);%����ֱ��������ھ�̬��ͷ����
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            end
                        end
                    end
                end
            else  %�����ڴ�ͷ         
                if(S3(i).f>S3(S3(i).NB(j)).f)
                     S3(i).Cid=S3(i).NB(j);%ѡ����ɼ�����������ھӽڵ�����
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                     S3(i).Cid='sink';%����ֱ��������ھ�̬��ͷ����
                end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%ѡ�����ߵȼ�
for i=1:1:n
    u=[];
    if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'����Ϊ����ڵ�
        for j=1:1:length(S3(i).S)
            u(j)=S3(S3(i).S(j)).M;
        end
        if(length(find(u==0))>0)
            if((S3(Sid(S3(i).C)).f>fo  && Sid(S3(i).C)~=0) || LB(i)==0)%���ڶ���̬��ͷ�����ӷ�ʽ
                S3(i).M=2;
            else
                S3(i).M=1;
            end
        end
    end
end
%���Ѵ������߽ڵ㣨��һ��ʼ�ھӽڵ��ǻ�Ծ�ģ����������ˣ�
for i=1:1:n
    u=[];
    if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'����Ϊ����ڵ�
        for j=1:1:length(S3(i).S)
            u(j)=S3(S3(i).S(j)).M;
        end
        if(length(find(u==0))==0)
            S3(i).M=0; 
        end
    end
end

for r=1:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r
    %�����ڵ���ģ��
dead3=0;
for i=1:1:n
    %������������ڵ�
    if (S3(i).E<=0)
        dead3=dead3+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead3==1)
           if(flag_first_dead3==0)
              first_dead3=r;
              flag_first_dead3=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead3==0.05*n)
           if(flag_teenth_dead3==0)%ÿ����5���ڵ�
              teenth_dead3=r;
              flag_teenth_dead3=1;
           end
        end
        if(dead3==n)
           if(flag_all_dead3==0)
              all_dead3=r;
              flag_all_dead3=1;
           end
        end
    end
    if S3(i).E>0
        S3(i).type='N';
    end
end
STATISTICS.DEAD3(r+1)=dead3;%�����ڵ�����
STATISTICS.ALLIVE3(r+1)=allive3-dead3;%���ڵ���
m=allive3-dead3;

%%%%%%%%%%%%%%%%%%%%%%%%%
%��վλ��ɸѡ
LB1=zeros(n,n);
SumLB1=zeros(1,n);
for j=1:1:n
    for i=1:1:n
       if (S3(i).E>0 && j~=i)
           d7=sqrt((S3(j).xd-S3(i).xd)^2 + (S5(j).yd-S5(i).yd)^2);
           LB(j,i)=d7;
       end
    end
end
SumLB1=sum(LB1,2);
[st7,BS]=min(SumLB1);

sink.x=Sink(r+1).x;
sink.y=Sink(r+1).y;

% sink.x=S3(BS).xd;
% sink.y=S3(BS).yd;

for i=1:1:n
    S3(i).D=sqrt((sink.x-S3(i).xd)^2 + (sink.y-S3(i).yd)^2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ھӽڵ��
for i=1:1:n
    if(S3(i).E>0)
        Contigupus=zeros(1,n);%�洢�ھӽڵ���
        Similar=zeros(1,n);
        A=zeros(1,n);
        for j=1:1:n
            d=sqrt((S3(i).xd-(S3(j).xd))^2 + (S3(i).yd-(S3(j).yd))^2);
            if(d<=2*L && d~=0 && S3(j).E>0)%�ڵ��������ڸ�֪�뾶��Ϊ�ھӽڵ�
                Contigupus(j)=j;   
            end
            if(d<L && d~=0 && S3(j).E>0)%��ѡ���ƽڵ�
                Similar(j)=j;   
            end
            if(S3(j).D>do && d<=do && d~=0 && S3(j).E>0 && S3(i).C==S3(j).C)
                A(j)=j;
            end
        end
        if (sum(find(Contigupus~=0))~=0)%�жϾ����Ƿ�Ϊ��
             S3(i).NB=find(Contigupus~=0);
        else
             S3(i).NB='null';
        end
        if (sum(find(Similar~=0))~=0)%�жϾ����Ƿ�Ϊ��
             S3(i).S=find(Similar~=0);
        else
             S3(i).S='null';
        end
        if (sum(find(A~=0))~=0)%�жϾ����Ƿ�Ϊ��
            S3(i).a=length(find(A~=0))+1;
        else
            S3(i).a=1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%�����ڵ����ڵ���75ʱ
for i=1:1:n
    if((75<m) && (m<=100))
         if (S3(i).yd>=S3(i).xd && S3(i).yd>(-S3(i).xd))
             S3(i).C=1;
        elseif (S3(i).yd<S3(i).xd && S3(i).yd>=(-S3(i).xd))
             S3(i).C=2;
        elseif (S3(i).yd<=S3(i).xd && S3(i).yd<(-S3(i).xd))
             S3(i).C=3;
        else
             S3(i).C=4;
        end
    elseif((50<m) && (m<=75))
        if (S3(i).yd>=0 && S3(i).yd>(-S3(i).xd))
             S3(i).C=1;
        elseif (S3(i).yd<=S3(i).xd && S3(i).yd<0)
             S3(i).C=2;
        else
             S3(i).C=3;
        end
    elseif((25<m) && (m<=50))
         if (S3(i).yd>=(-S3(i).xd))
             S3(i).C=1;
         else 
             S3(i).C=2;
         end
    else
             S3(i).C=1;
    end
end
  %%%%%%%%%%%%%%%%%%%%%
   %����С�㼶����
e=zeros(4,B);%4��ʾ������
a=zeros(4,B+1);
eav=zeros(4,B+1);
for j=1:1:n
   e(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)+S3(j).E;%����С�㼶������
   a(S3(j).C,S3(j).f)=a(S3(j).C,S3(j).f)+1;%���������С�㼶�Ľڵ���
end
for j=1:1:n
   eav(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)./a(S3(j).C,S3(j).f);
   a(S3(j).C,B+1)= a(S3(j).C,B+1)+1;%���������Ľڵ����
end
Et=0;
for i=1:1:n
    if(S3(i).E>0)
        Et=Et+S3(i).E;%��ʣ������
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ea=Et/(allive3-dead3);%EaΪʣ����ƽ������
SCHCount=ceil(0.04*m);
N=zeros(SCHCount,n);%���澲̬��ͷɸѡ����
R=zeros(1,n);%���涯̬��ͷɸѡ����
JE=zeros(1,SCHCount);%�����ж������Ƿ���ڶ�̬��ͷ
countCHs3=0;%�����ͷ����ɸѡ��������ͷ��������
if(m>=10)
for i=1:1:n
   N(S3(i).C,i)=(eav(S3(i).C,S3(i).f)*(S3(i).E-Ea))*(S3(i).a/n)/(4*Eo-S3(i).E)*(1-(S3(i).D/(1.42*B)));
%  if Ea>0
%      p(i)=P*S3(i).a/n;
%      if(S3(i).E>0)    
%            %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
%         N(S3(i).C,i)=(p(i)/(1-p(i)*mod(r,round(1/p(i))))*(S3(i).E/(3*eav(S3(i).C,S3(i).f))));%p(i)<1,�����
%      end  
%   end 
end
[st,Sid]=max(N,[],2);%��ȡ��̬��ͷɸѡ�����е����ֵ����ID
for j=1:1:SCHCount
if (st(j)~=0 && a(S3(Sid(j)).C,B+1)>10) %��������������������̬��ͷ
   S3(Sid(j)).type='C';
   S3(Sid(j)).Cid='sink';%û�ж�̬��ͷʱ����̬��ͷ���ڵ�sink����
   countCHs3=countCHs3+1;
   if( S3(Sid(j)).f>fo )
       for z=1:1:n 
           if(S3(z).D<S3(Sid(j)).D && S3(z).C==S3(Sid(j)).C)%ͬһ����͵���ڵ�������ھ�̬��ͷʱ���㶯̬�ڵ�Ӱ������
               R(z)=(4*(S3(Sid(j)).f-S3(z).f+1)*S3(z).f/((S3(Sid(j)).f+1)*(S3(Sid(j)).f+1)))*(1-S3(z).f/S3(Sid(j)).f)*(eav(S3(z).C,S3(z).f)*(S3(j).E-Ea))/(3*Eo-S3(z).E);
           end
       end
   [dy,Did]=max(R,[],2); %��ȡ��̬��ͷɸѡ�����е����ֵ����ID  
   S3(Did).type='D';
   S3(Sid(j)).Cid=Did;%�ж�̬��ͷ������Ϊ�붯̬��ͷ����
   S3(Did).Cid='sink';%��̬��ͷ���ڵ�����
   JE(S3(Did).C)=Did;
   countCHs3=countCHs3+1;
   end
end
end
end
STATISTICS.COUNTCHS3(r+1)=countCHs3;
%%%%%%%%%%%%%%%%%%%%%%
%�ڵ�����
LB=zeros(1,n);%������Ϣ������
for i=1:1:n
    if ( S3(i).type=='N' && S3(i).E>0 && S3(i).D<=do)
        S3(i).Cid='sink';
    elseif(S3(i).type=='N' && S3(i).E>0 && S3(i).D>do)
            if(JE(S3(i).C)~=0  && Sid(S3(i).C)~=0)%���ڶ�����̬��ͷ
                d1=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                d2=sqrt((S3(i).xd-S3(JE(S3(i).C)).xd)^2+(S3(i).yd-S3(JE(S3(i).C)).yd)^2);
                if(d1<d2)
                    S3(i).Cid=Sid(S3(i).C);%����ڵ������͵�����ڵ�ID
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                    S3(i).Cid=JE(S3(i).C);
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                end
            elseif(Sid(S3(i).C)~=0 && JE(S3(i).C)==0)%�����ڶ�̬��ͷ�����ھ�̬��ͷ
                d3=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                if(sum(S3(i).NB)==443 || d3<=do)%���û���ھӽڵ㣬��ֱ���뾲̬��ͷ����,sum('null')==443
                     S3(i).Cid=Sid(S3(i).C);
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                elseif(sum(S3(i).NB)~=443 && d3>do)
                    b=[];
                    for z=1:1:length(S3(i).NB)
                        b(z)=S3(S3(i).NB(z)).E;
                    end
                    for j=1:1:length(S3(i).NB)
                        if(S3(S3(i).NB(j)).type~='N')%���ھӽڵ���������ڵ�ʱ,��������ڵ�
                            S3(i).Cid=S3(i).NB(j);
                            LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            break;
                        else
                            if(S3(i).f>S3(S3(i).NB(j)).f && S3(S3(i).NB(j)).E==max(b))
                                 S3(i).Cid=S3(i).NB(j);%ѡ����㼶�������������������ھӽڵ�����
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                                 break;
                            else
                                 S3(i).Cid=Sid(S3(i).C);%����ֱ��������ھ�̬��ͷ����
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            end
                        end
                    end
                end
            else  %�����ڴ�ͷ         
                if(S3(i).f>S3(S3(i).NB(j)).f)
                     S3(i).Cid=S3(i).NB(j);%ѡ����ɼ�����������ھӽڵ�����
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                     S3(i).Cid='sink';%����ֱ��������ھ�̬��ͷ����
                end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%ѡ�����ߵȼ�
for i=1:1:n
    if(S3(i).E>0)
        u=[];
        if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'����Ϊ����ڵ�
            for j=1:1:length(S3(i).S)
                u(j)=S3(S3(i).S(j)).M;
            end
            if(length(find(u==0))>0)
                if((S3(Sid(S3(i).C)).f>fo  && Sid(S3(i).C)~=0) || LB(i)==0)%���ڶ���̬��ͷ�����ӷ�ʽ
                    S3(i).M=2;
                else
                    S3(i).M=1;
                end
            end
        end
    else
        S3(i).M=2;
    end
end
%���Ѵ������߽ڵ㣨��һ��ʼ�ھӽڵ��ǻ�Ծ�ģ����������ˣ�
for i=1:1:n
     if(S3(i).E>0)
        u=[];
        if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'����Ϊ����ڵ�
            for j=1:1:length(S3(i).S)
                u(j)=S3(S3(i).S(j)).M;
            end
            if(length(find(u==0))==0)
                S3(i).M=0; 
            end
        end
        if(sum(S3(i).Cid)~=437 && S3(i).M~=0)%��ȥ���߽ڵ��������
            LB(S3(i).Cid)=LB(S3(i).Cid)-1;
        end
     end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Ϊ��ͷ��������
% F=zeros(1,n);
% for i=1:1:n
%     a=1;
%     g=Inf;
%     Z0 = zeros(1,n);
%     if(sum(S3(i).S)~=443)
%         for j=1:1:length(S3(i).S)
%             if (S3(S3(i).S(j)).M==2 && S3(S3(i).S(j)).E > 0)
%                 a=a+1;
%                 if(g>S3(S3(i).S(j)).E)
%                     g = S3(S3(i).S(j)).E;
%                 end
%             end
%         end
%        F(i) = a/g;
%     end
% end
% [w,x]=max(F);
% if(sum(S3(i).S)~=443)
%     for i=1:1:length(S3(w).S)
%          if (S3(S3(w).S(i)).M==2 && S3(S3(w).S(i)).E > 0)
%              S3(S3(w).S(i)).E=0.01*S3(S3(w).S(i)).D+S3(S3(w).S(i)).E;
%          end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ڵ����ϼ��ڵ�Cid�ľ���
DCid3 = zeros(1,B);
for i=1:1:n
    %if (S3(i).M~=2)
        if(sum(S3(i).Cid)~=437)
            DCid3(1,i)=sqrt((S3(i).xd-S3(S3(i).Cid).xd)^2+(S3(i).yd-S3(S3(i).Cid).yd)^2); 
        else
            DCid3(1,i)=S3(i).D;
        end
        AS3(1,q)=AS3(1,q)+1;
    %end
end
for i=1:1:n
    if (DCid3(i)~=0 && DCid3(i)>=do)
        if(((PL+LB(i))*rand(1,1)+(DCid3(i)/do)-1)>PL)
            Loss3(1,q)=Loss3(1,q)+1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݴ����������ݽ�����Ϊ100%
for i=1:1:n
    if(S3(i).E>0 && S3(i).M~=2)
       data3=data3+1; 
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %�ܺ�ģ��
%�ڵ�ÿ�ַ���һ�����ݣ�ÿ��Ϊ4000bit
for i=1:1:n
    S3(i).EB=S3(i).E;
end
for i=1:1:n
   if (S3(i).M==0)
        if (sum(S3(i).Cid)~=437)
            d1=sqrt((S3(i).xd-S3(S3(i).Cid).xd)^2+(S3(i).yd-S3(S3(i).Cid).yd)^2);
        if (d1>do)
            S3(i).E=S3(i).E-((EPC+ETX)*4000+(EDA+ERX)*4000*LB(i)+Emp*4000*( d1 * d1 * d1 * d1)); 
        end
        if (d1<=do)
            S3(i).E=S3(i).E-((EPC+ETX)*4000+(EDA+ERX)*4000*LB(i)+Efs*4000*( d1 * d1)); 
        end
        else
          d1=sqrt(S3(i).xd^2+S3(i).yd^2);
          if (d1>do)
            S3(i).E=S3(i).E-((EPC+ETX)*4000+(EDA+ERX)*4000*LB(i)+Emp*4000*( d1 * d1 * d1 * d1)); 
        end
        if (d1<=do)
            S3(i).E=S3(i).E-((EPC+ETX)*4000+(EDA+ERX)*4000*LB(i)+Efs*4000*( d1 * d1)); 
        end
        end
   elseif(S3(i).M==1)
        if (sum(S3(i).Cid)~=437)
        d1=sqrt((S3(i).xd-S3(S3(i).Cid).xd)^2+(S3(i).yd-S3(S3(i).Cid).yd)^2);
        if (d1>do)
            S3(i).E=S3(i).E-((ETX)*4000+(EDA+ERX)*4000*LB(i)+Emp*4000*( d1 * d1 * d1 * d1)); 
        end
        if (d1<=do)
            S3(i).E=S3(i).E-((ETX)*4000+(EDA+ERX)*4000*LB(i)+Efs*4000*( d1 * d1)); 
        end
        else
          d1=sqrt(S3(i).xd^2+S3(i).yd^2);
          if (d1>do)
            S3(i).E=S3(i).E-((ETX)*4000+(EDA+ERX)*4000*LB(i)+Emp*4000*( d1 * d1 * d1 * d1)); 
        end
        if (d1<=do)
            S3(i).E=S3(i).E-((ETX)*4000+(EDA+ERX)*4000*LB(i)+Efs*4000*( d1 * d1)); 
        end
        end
   end
end
if(countCHs3>0)%�д�ͷʱ����
    EB3=[];
    Eb3=1;
    Cb3=0;
    for i=1:1:n
        if (S3(i).type=='C' || S3(i).type=='D')%���д�ͷ������������ʱ
            if (S3(i).EB>0 && S3(i).E>0)
                EB3(Eb3)=(S3(i).EB-S3(i).E)/S3(i).EB*100;
                if (EB3(Eb3)>Cb3)
                    Cb3=EB3(Eb3);
                end
                Eb3=Eb3+1;
            end
        end
    end
    if(length(EB3)>=1)%���д�ͷ��������������ʱ
        Ave3=sum(EB3)/length(EB3);
        STATISTICS.Ave3(r+1)=Ave3;
        STATISTICS.Cb3(r+1)=Cb3;
    else
        STATISTICS.Ave3(r+1)=0;
        STATISTICS.Cb3(r+1)=0;
    end
else
    STATISTICS.Ave3(r+1)=0;
    STATISTICS.Cb3(r+1)=0;
end
% %%%%%%%%%%%
% %����������
for i=1:1:n
    if(S3(i).E>0)
        S3(i).M=0;
    else
        S3(i).M=2;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leach�㷨
for r=1:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r
  %ÿ��һ����ת����(������Ϊ20��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣�
  % �ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ�㣬��50�����ڣ�ÿ����20��
  if(mod(r,round(1/P))==0)%round��������ȡ������, mod(x,y)Ϊ���ຯ����yΪ�����������������ͬ
    for i=1:1:n
        S4(i).G=0;
    end
  end
  
  Et4=0;
  for i=1:1:n
      if(S4(i).E>0)
        Et4=S4(i).E+Et4;%El3Ϊʣ��������
      end
  end
%(2)�����ڵ���ģ��
dead4=0;
for i=1:1:n
    %������������ڵ�
    if (S4(i).E<=0)
        dead4=dead4+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead4==1)
           if(flag_first_dead4==0)
              first_dead4=r;
              flag_first_dead4=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead4==0.1*n)
           if(flag_teenth_dead4==0)%�����һ�����ڣ�20��
              teenth_dead4=r;
              flag_teenth_dead4=1;
           end
        end
        if(dead4==n)
           if(flag_all_dead4==0)
              all_dead4=r;
              flag_all_dead4=1;
           end
        end
    end
    if S4(i).E>0
        S4(i).type='N';
    end
end
STATISTICS.DEAD4(r+1)=dead4;%�����ڵ�����
STATISTICS.ALLIVE4(r+1)=allive4-dead4;%���ڵ���
m4=allive4-dead4;
Ea4=Et4/(allive4-dead4);%EtΪ��������Et*(1-r/rmax)Ϊÿһ��ʣ�����������EaΪʣ����ƽ������
%(4)��ͷѡ��ģ��
countCHs4=0;
cluster4=1;
C4=[];
for i=1:1:n
%  if Ea4>0
%  p(i)=P*n*S4(i).E*S4(i).E/(Et*Ea4);% E3(i)= S3(i).E
 if(S4(i).E>0 && Ea4>0)
   temp_rand=rand;     
   if ( (S4(i).G)<=0)  
       %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
        if(temp_rand<= (P/(1-P*mod(r,round(1/P)))))%p(i)<1,�����
            %p(i)/(1-p(i)*mod(r,round(1/p(i)))))Ϊ��ֵ���㹫ʽ
            countCHs4=countCHs4+1;
            S4(i).type='C';
            S4(i).Cid='sink';
            S4(i).G=round(1/P)-1;%p(i)ԽС��S3(i).GԽ�󣬴���0ʱ���ܾ�ѡ��ͷ
            C4(cluster4)=i;
            cluster4=cluster4+1;
        end      
    end
    % S3(i).G=S3(i).G-1;     
  end 
%  end
end
STATISTICS.COUNTCHS4(r+1)=countCHs4;%��ͷ����

%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
LB4=zeros(1,n);
for i=1:1:n
   if ( S4(i).type=='N' && S4(i).E>0 )
       min_dis=Inf;%iInf��ʾ�������
       for c=1:1:cluster4-1
           temp=min(min_dis,sqrt((S4(i).xd-S4(C4(c)).xd)^2 + (S4(i).yd-S4(C4(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S4(i).Cid=C4(c);%������ͨ�����ӵĴ�ͷ�ı��
               LB4(C4(c))=LB4(C4(c))+1;%������ڽڵ����
           end
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ڵ����ϼ��ڵ�Cid�ľ���
DCid4 = zeros(1,B);
for i=1:1:n
    if (S4(i).E>0)
        if(sum(S4(i).Cid)~=437)
            DCid4(1,i)=sqrt((S4(i).xd-S4(S4(i).Cid).xd)^2+(S4(i).yd-S4(S4(i).Cid).yd)^2); 
        else
            DCid4(1,i)=S4(i).D;
        end
    end
end
for i=1:1:n
    if (DCid4(i)~=0 && DCid4(i)>=do)
        if(((PL+LB4(i))*rand(1,1)+(DCid4(i)/do)-1)>PL)
            Loss4(1,q)=Loss4(1,q)+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݴ����������ݽ�����Ϊ100%
for i=1:1:n
    if(S4(i).E>0)
       data4=data4+1; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%�ܺ�ģ��
for i=1:1:n
    S4(i).EB=S4(i).E;
end
 for i=1:1:n
     if(S4(i).E>0)
         if(sum(S4(i).Cid)~=437 )
             d=sqrt((S4(i).xd-S4(S4(i).Cid).xd)^2 + (S4(i).yd-S4(S4(i).Cid).yd)^2);
            if (d>do)
                    S4(i).E=S4(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB4(i) + Emp*4000*( d * d * d * d));
                end
                if (d<=do)
                    S4(i).E=S4(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB4(i) + Efs*4000*( d * d)); 
                end
         else
              d=sqrt((S4(i).xd)^2 + (S4(i).yd)^2);
              if (d>do)
                    S4(i).E=S4(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB4(i) + Emp*4000*( d * d * d * d)); 
                end
                if (d<=do)
                    S4(i).E=S4(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB4(i) + Efs*4000*( d * d)); 
                end
         end
     end
 end
if(countCHs4>0)%�д�ͷʱ����
    EB4=[];
    Eb4=1;
    Cb4=0;
    for i=1:1:n
        if (S4(i).type=='C')%���д�ͷ������������ʱ
            if (S4(i).EB>0 && S4(i).E>0)
                EB4(Eb4)=(S4(i).EB-S4(i).E)/S4(i).EB*100;
                if (EB4(Eb4)>Cb4)
                    Cb4=EB4(Eb4);
                end
                Eb4=Eb4+1;
            end
        end
    end
    if(length(EB4)>=1)%���д�ͷ��������������ʱ
        Ave4=sum(EB4)/length(EB4);
        STATISTICS.Ave4(r+1)=Ave4;
        STATISTICS.Cb4(r+1)=Cb4;
    else
        STATISTICS.Ave4(r+1)=0;
        STATISTICS.Cb4(r+1)=0;
    end
else
    STATISTICS.Ave4(r+1)=0;
    STATISTICS.Cb4(r+1)=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEACH-N�㷨
for r=1:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r
  %ÿ��һ����ת����(������Ϊ20��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣�
  % �ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ�㣬��50�����ڣ�ÿ����20��
  if(mod(r,round(1/P))==0)%round��������ȡ������, mod(x,y)Ϊ���ຯ����yΪ�����������������ͬ
    for i=1:1:n
        S5(i).G=0;
    end
  end 
  Et5=0;
  for i=1:1:n
      if(S5(i).E>0)
        Et5=S5(i).E+Et5;%El3Ϊʣ��������
      end
  end
%(2)�����ڵ���ģ��
dead5=0;
for i=1:1:n
    %������������ڵ�
    if (S5(i).E<=0)
        dead5=dead5+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead5==1)
           if(flag_first_dead5==0)
              first_dead5=r;
              flag_first_dead5=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead5==0.1*n)
           if(flag_teenth_dead5==0)%�����һ�����ڣ�20��
              teenth_dead5=r;
              flag_teenth_dead5=1;
           end
        end
        if(dead5==n)
           if(flag_all_dead5==0)
              all_dead5=r;
              flag_all_dead5=1;
           end
        end
    end
    if S5(i).E>0
        S5(i).type='N';
    end
end
STATISTICS.DEAD5(r+1)=dead5;%�����ڵ�����
STATISTICS.ALLIVE5(r+1)=allive5-dead5;%���ڵ���
m5=allive5-dead5;
Ea5=Et5/(allive5-dead5);%EtΪ��������Et*(1-r/rmax)Ϊÿһ��ʣ�����������EaΪʣ����ƽ������

%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ھӽڵ��
for i=1:1:n
    if(S5(i).E>0)
        A=zeros(1,n);%�洢�ھӽڵ���
        for j=1:1:n
            d=sqrt((S5(i).xd-(S5(j).xd))^2 + (S5(i).yd-(S5(j).yd))^2);
            if(d<=do && d~=0 && S5(j).E>0)%�ڵ��������ڸ�֪�뾶��Ϊ�ھӽڵ�
                A(j)=j;   
            end
        end
        if (sum(find(A~=0))~=0)%�жϾ����Ƿ�Ϊ��
             S5(i).a=length(find(A~=0));
        else
             S5(i).a=0;
        end
    end
end

%(4)��ͷѡ��ģ��
countCHs5=0;
cluster5=1;
C5=[];
for i=1:1:n
 if Ea5>0
 p(i)=P*S5(i).a/n;
 if(S5(i).E>0)
   temp_rand=rand;     
   if ( (S5(i).G)<=0)  
       %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))*(S5(i).E/(3*Eo))))%p(i)<1,�����
            %p(i)/(1-p(i)*mod(r,round(1/p(i)))))Ϊ��ֵ���㹫ʽ
            countCHs5=countCHs5+1;
            S5(i).type='C';
            S5(i).Cid='sink';
            S5(i).G=round(1/p(i))-1;%p(i)ԽС��S3(i).GԽ�󣬴���0ʱ���ܾ�ѡ��ͷ
            C5(cluster5)=i;
            cluster5=cluster5+1;
        end      
    end
    % S3(i).G=S3(i).G-1;     
  end 
 end
end
STATISTICS.COUNTCHS5(r+1)=countCHs5;%��ͷ����

%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
LB5=zeros(1,n);
for i=1:1:n
   if ( S5(i).type=='N' && S5(i).E>0 )
       min_dis=Inf;%iInf��ʾ�������
       for c=1:1:cluster5-1
           temp=min(min_dis,sqrt((S5(i).xd-S5(C5(c)).xd)^2 + (S5(i).yd-S5(C5(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S5(i).Cid=C5(c);%������ͨ�����ӵĴ�ͷ�ı��
               LB5(S5(i).Cid)=LB5(S5(i).Cid)+1;%������ڽڵ����
           end
       end
%        if(min_dis<S5(i).D)
%            LB5(S5(i).Cid)=LB5(S5(i).Cid)+1;%������ڽڵ����
%        else
%            S5(i).Cid='sink';
%        end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ڵ����ϼ��ڵ�Cid�ľ���
DCid5 = zeros(1,B);
for i=1:1:n
    if (S5(i).E>0)
        if(sum(S5(i).Cid)~=437)
            DCid5(1,i)=sqrt((S5(i).xd-S5(S5(i).Cid).xd)^2+(S5(i).yd-S5(S5(i).Cid).yd)^2); 
        else
            DCid5(1,i)=S5(i).D;
        end
    end
end
for i=1:1:n
    if (DCid5(i)~=0 && DCid5(i)>=do)
        if(((PL+LB5(i))*rand(1,1)+(DCid5(i)/do)-1)>PL)
            Loss5(1,q)=Loss5(1,q)+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݴ����������ݽ�����Ϊ100%
for i=1:1:n
    if(S5(i).E>0)
       data5=data5+1; 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�ܺ�ģ��
for i=1:1:n
    S5(i).EB=S5(i).E;
end
 for i=1:1:n
     if(S5(i).E>0)
         if(sum(S5(i).Cid)~=437 )
             d=sqrt((S5(i).xd-S5(S5(i).Cid).xd)^2 + (S5(i).yd-S5(S5(i).Cid).yd)^2);
            if (d>do)
                    S5(i).E=S5(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i) + Emp*4000*( d * d * d * d));
                end
                if (d<=do)
                    S5(i).E=S5(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB5(i) + Efs*4000*( d * d)); 
                end
         else
              d=sqrt((S5(i).xd)^2 + (S5(i).yd)^2);
              if (d>do)
                    S5(i).E=S5(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB5(i) + Emp*4000*( d * d * d * d)); 
                end
                if (d<=do)
                    S5(i).E=S5(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB5(i) + Efs*4000*( d * d)); 
                end
         end
     end
 end
if(countCHs5>0)%�д�ͷʱ����
    EB5=[];
    Eb5=1;
    Cb5=0;
    for i=1:1:n
        if (S5(i).type=='C')%���д�ͷ������������ʱ
            if (S5(i).EB>0 && S5(i).E>0)
                EB5(Eb5)=(S5(i).EB-S5(i).E)/S5(i).EB*100;
                if (EB5(Eb5)>Cb5)
                    Cb5=EB5(Eb5);
                end
                Eb5=Eb5+1;
            end
        end
    end
    if(length(EB5)>=1)%���д�ͷ��������������ʱ
        Ave5=sum(EB5)/length(EB5);
        STATISTICS.Ave5(r+1)=Ave5;
        STATISTICS.Cb5(r+1)=Cb5;
    else
        STATISTICS.Ave5(r+1)=0;
        STATISTICS.Cb5(r+1)=0;
    end
else
    STATISTICS.Ave5(r+1)=0;
    STATISTICS.Cb5(r+1)=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEACH-C�㷨
for r=1:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r
  %ÿ��һ����ת����(������Ϊ20��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣�
  % �ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ�㣬��50�����ڣ�ÿ����20��
  if(mod(r,round(1/P))==0)%round��������ȡ������, mod(x,y)Ϊ���ຯ����yΪ�����������������ͬ
    for i=1:1:n
        S6(i).G=0;
    end
  end
  
  Et6=0;
for i=1:100
   if(S6(i).E>0)
        Et6=S6(i).E+Et6;%El3Ϊʣ��������
   end
end
%(2)�����ڵ���ģ��
dead6=0;
for i=1:1:n
    %������������ڵ�
    if (S6(i).E<=0)
        dead6=dead6+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead6==1)
           if(flag_first_dead6==0)
              first_dead6=r;
              flag_first_dead6=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead6==0.1*n)
           if(flag_teenth_dead6==0)%�����һ�����ڣ�20��
              teenth_dead6=r;
              flag_teenth_dead6=1;
           end
        end
        if(dead6==n)
           if(flag_all_dead6==0)
              all_dead6=r;
              flag_all_dead6=1;
           end
        end
    end
    if S6(i).E>0
        S6(i).type='N';
    end
end
STATISTICS.DEAD6(r+1)=dead6;%�����ڵ�����
STATISTICS.ALLIVE6(r+1)=allive6-dead6;%���ڵ���
m6=allive6-dead6;
Ea6=Et6/(allive6-dead6);%EtΪ��������Et*(1-r/rmax)Ϊÿһ��ʣ�����������EaΪʣ����ƽ������

%%%%%%%%%%%%%%%%%%%%%%%%
%�����ڵ����ڵ���75ʱ
for i=1:1:n
    if((75<m6) && (m6<=100))
         if (S6(i).yd>=S6(i).xd && S6(i).yd>(-S6(i).xd))
             S6(i).C=1;
        elseif (S6(i).yd<S6(i).xd && S6(i).yd>=(-S6(i).xd))
             S6(i).C=2;
        elseif (S6(i).yd<=S6(i).xd && S6(i).yd<(-S6(i).xd))
             S6(i).C=3;
        else
             S6(i).C=4;
        end
    elseif((50<m6) && (m6<=75))
        if (S6(i).yd>=0 && S6(i).yd>(-S6(i).xd))
             S6(i).C=1;
        elseif (S6(i).yd<=S6(i).xd && S6(i).yd<0)
             S6(i).C=2;
        else
             S6(i).C=3;
        end
    elseif((25<m6) && (m6<=50))
         if (S6(i).yd>=(-S6(i).xd))
             S6(i).C=1;
         else 
             S6(i).C=2;
         end
    else
             S6(i).C=1;
    end
end

%(4)��ͷѡ��ģ��
SCHCount6=ceil(0.04*m6);
countCHs6=0;
cluster6=1;
C6=[];
D=zeros(n,n);
DL=[];
for i=1:1:n
    if(Ea6>0 && S6(i).E>Ea)
        for z=1:1:n
            if(S6(i).C==S6(z).C && i~=z)
                D(i,z)=sqrt(sqrt((S6(i).xd-(S6(z).xd))^2 + (S6(i).yd-(S6(z).yd))^2));
            end
        end
    end
end
DL=sum(D,2);
for i=1:1:n
    if(DL(i)==0)
       DL(i)=Inf;
    end
end
for i=1:1:SCHCount6
    [j,k]=min(DL);
    S6(k).type='C';
    S6(k).Cid='sink';
    C6(cluster6)=k;
    cluster6=cluster6+1;
    countCHs6=countCHs6+1;
   for z=1:1:n
       if(S6(k).C==S6(z).C)
           DL(z)=Inf;
       end
   end
end   
STATISTICS.COUNTCHS6(r+1)=countCHs6;%��ͷ����
%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
LB6=zeros(1,n);
for i=1:1:n
   if ( S6(i).type=='N' && S6(i).E>0 )
       min_dis=Inf;%iInf��ʾ�������
       for c=1:1:cluster6-1
           temp=min(min_dis,sqrt((S6(i).xd-S6(C6(c)).xd)^2 + (S6(i).yd-S6(C6(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S6(i).Cid=C6(c);%������ͨ�����ӵĴ�ͷ�ı��
               LB6(C6(c))=LB6(C6(c))+1;%������ڽڵ����
           end
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ڵ����ϼ��ڵ�Cid�ľ���
DCid6 = zeros(1,B);
for i=1:1:n
    if (S6(i).E>0)
        if(sum(S6(i).Cid)~=437)
            DCid6(1,i)=sqrt((S6(i).xd-S6(S6(i).Cid).xd)^2+(S6(i).yd-S6(S6(i).Cid).yd)^2); 
        else
            DCid6(1,i)=S6(i).D;
        end
    end
end
for i=1:1:n
    if (DCid6(i)~=0 && DCid6(i)>=do)
        if(((PL+LB6(i))*rand(1,1)+(DCid6(i)/do)-1)>PL)
            Loss6(1,q)=Loss6(1,q)+1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݴ����������ݽ�����Ϊ100%
for i=1:1:n
    if(S6(i).E>0)
       data6=data6+1; 
    end
end

for i=1:1:n
    S6(i).EB=S6(i).E;
end
for i=1:1:n
     if(S6(i).E>0)
         if(sum(S6(i).Cid)~=437 )
             d=sqrt((S6(i).xd-S6(S6(i).Cid).xd)^2 + (S6(i).yd-S6(S6(i).Cid).yd)^2);
            if (d>do)
                S6(i).E=S6(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB6(i) + Emp*4000*( d * d * d * d));
            end
            if (d<=do)
                S6(i).E=S6(i).E-( (EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB6(i) + Efs*4000*( d * d)); 
            end
         else
              d=sqrt((S6(i).xd)^2 + (S6(i).yd)^2);
            if (d>do)
                S6(i).E=S6(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB6(i) + Emp*4000*( d * d * d * d)); 
            end
            if (d<=do)
                S6(i).E=S6(i).E-( ETX*(4000)+(ERX+EDA)*(4000)*LB6(i) + Efs*4000*( d * d)); 
            end
         end
     end
end
if(countCHs6>0)%�д�ͷʱ����
    EB6=[];
    Eb6=1;
    Cb6=0;
    for i=1:1:n
        if (S6(i).type=='C')%���д�ͷ������������ʱ
            if (S6(i).EB>0 && S6(i).E>0)
                EB6(Eb6)=(S6(i).EB-S6(i).E)/S6(i).EB*100;
                if (EB6(Eb6)>Cb6)
                    Cb6=EB6(Eb6);
                end
                Eb6=Eb6+1;
            end
        end
    end
    if(length(EB6)>=1)%���д�ͷ��������������ʱ
        Ave6=sum(EB6)/length(EB6);
        STATISTICS.Ave6(r+1)=Ave6;
        STATISTICS.Cb6(r+1)=Cb6;
    else
        STATISTICS.Ave6(r+1)=0;
        STATISTICS.Cb6(r+1)=0;
    end
else
    STATISTICS.Ave6(r+1)=0;
    STATISTICS.Cb6(r+1)=0;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEACH-R�㷨
for r=1:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r 
BS.x=(rand(1,1)*2-1)*xm;
BS.y=(rand(1,1)*2-1)*ym;
for i=1:1:n
    S7(i).D=sqrt((S7(i).xd-BS.x)^2 + (S7(i).yd-BS.y)^2);
end
  %ÿ��һ����ת����(������Ϊ20��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣�
  % �ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ�㣬��50�����ڣ�ÿ����20��
  if(mod(r,round(1/P))==0)%round��������ȡ������, mod(x,y)Ϊ���ຯ����yΪ�����������������ͬ
    for i=1:1:n
        S7(i).G=0;
    end
  end
  Et7=0;
  for i=1:1:n
      if(S7(i).E>0)
        Et7=S7(i).E+Et7;%El3Ϊʣ��������
      end
  end
%(2)�����ڵ���ģ��
dead7=0;
for i=1:1:n
    %������������ڵ�
    if (S7(i).E<=0)
        dead7=dead7+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead7==1)
           if(flag_first_dead7==0)
              first_dead7=r;
              flag_first_dead7=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead7==0.1*n)
           if(flag_teenth_dead7==0)%�����һ�����ڣ�20��
              teenth_dead7=r;
              flag_teenth_dead7=1;
           end
        end
        if(dead7==n)
           if(flag_all_dead7==0)
              all_dead7=r;
              flag_all_dead7=1;
           end
        end
    end
    if S7(i).E>0
        S7(i).type='N';
    end
end
STATISTICS.DEAD7(r+1)=dead7;%�����ڵ�����
STATISTICS.ALLIVE7(r+1)=allive7-dead7;%���ڵ���
m7=allive7-dead7;
Ea7=Et7/(allive7-dead7);%EtΪ��������Et*(1-r/rmax)Ϊÿһ��ʣ�����������EaΪʣ����ƽ������
%(4)��ͷѡ��ģ��
countCHs7=0;
cluster7=1;
C7=[];
Ymed=0;
if(r<=20)
    for i=1:1:n
     if(Ea7>0)
     if(S7(i).E>0)
       temp_rand=rand;     
       if ( (S7(i).G)<=0)  
           %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
            if(temp_rand<= (P/(1-P*mod(r,round(1/P)))))%p(i)<1,�����
                %p(i)/(1-p(i)*mod(r,round(1/p(i)))))Ϊ��ֵ���㹫ʽ
                countCHs7=countCHs7+1;
                S7(i).type='C';
                S7(i).Cid='BS';
                S7(i).G=round(1/P)-1;%p(i)ԽС��S7(i).GԽ�󣬴���0ʱ���ܾ�ѡ��ͷ
                C7(cluster7)=i;
                cluster7=cluster7+1;
            end      
        end
        % S3(i).G=S3(i).G-1;     
      end 
      Yn=0;
      for i=1:1:n
          if (sum(S7(i).NB)~=443)
              for j=1:1:length(S7(i).NB)
                  if (S7(i).E>S7(S7(i).NB(j)).E)
                        Yn=Yn+1;
                  else
                        Yn=Yn-1;
                  end
              end 
          end
          S7(i).Yn=Yn;
        end
     end
    end
else
   YN=[];
   for i=1:1:n
       YN(i)=S7(i).Yn;
   end
   Ymed=median(YN);
   for i=1:1:n
     if(S7(i).E>0 && Ea7>0 && S7(i).G<=0)
       temp_rand=rand;     
       if (Ymed~=0)  
           %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
            if (temp_rand <= (P/(1-P*mod(r,round(1/P))))*S7(i).Yn/(abs(Ymed)))%p(i)<1,�����
                countCHs7=countCHs7+1;
                S7(i).type='C';
                S7(i).Cid='BS';
                S7(i).G=round(1/P)-1;%p(i)ԽС��S7(i).GԽ�󣬴���0ʱ���ܾ�ѡ��ͷ
                C7(cluster7)=i;
                cluster7=cluster7+1;
            end 
       else
           if (temp_rand <= (P/(1-P*mod(r,round(1/P))))*S7(i).Yn)%p(i)<1,�����
                countCHs7=countCHs7+1;
                S7(i).type='C';
                S7(i).Cid='BS';
                S7(i).G=round(1/P)-1;
                C7(cluster7)=i;
                cluster7=cluster7+1;
           end  
       end  
    end
   Yn=0;
   for i=1:1:n
       if (sum(S7(i).NB)~=443)
           for j=1:1:length(S7(i).NB)
               if (S7(i).E>S7(S7(i).NB(j)).E)
                    Yn=Yn+1;
               else
                    Yn=Yn-1;
               end
           end 
       end
       S7(i).Yn=Yn;
   end
   end
end
STATISTICS.COUNTCHS7(r+1)=countCHs7;%��ͷ����

%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
LB7=zeros(1,n);
ds=Inf;
Dm=0;

for i=1:1:n
   if(S7(i).E>0)
   if ( S7(i).type=='N')
       min_dis=Inf;%iInf��ʾ�������
       for c=1:1:cluster7-1
           temp=min(min_dis,sqrt((S7(i).xd-S7(C7(c)).xd)^2 + (S7(i).yd-S7(C7(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S7(i).Cid=C7(c);%������ͨ�����ӵĴ�ͷ�ı��
               LB7(C7(c))=LB7(C7(c))+1;%������ڽڵ����
           end
       end
   end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:1:n
     if(S7(i).E>0)
     if (S7(i).type=='C')
         if(S7(i).D>do)
             for j=1:1:n %���������нڵ�ds������Сֵ
                 Dm=sqrt((S7(i).xd-S7(j).xd)^2 + (S7(i).yd-S7(j).yd)^2)+sqrt((S7(j).xd-BS.x)^2 + (S7(j).yd-BS.y)^2);
                 if(Dm<ds)
                     ds=Dm;
                     S7(i).Cid=j;
                     J=j;
                 end
             end
             if(S7(J).D<=do)
                 S7(J).Cid='BS';
                 LB7(S7(i).Cid)=LB7(S7(i).Cid)+1; 
             else
                 for j=1:1:n %��ȫ���ڵ�
                     if(J~=j)%S7(i).Cid~=j�ų��Լ�����û�ų��Լ����¼���ͷ
                         if(j~=i)%�ų��Լ����¼���ͷ
                         Dm=sqrt((S7(J).xd-S7(j).xd)^2 + (S7(J).xd-S7(j).yd)^2)+sqrt((S7(j).xd-BS.x)^2 + (S7(j).yd-BS.y)^2);
                            if(Dm<ds)
                                ds=Dm;
                                S7(J).Cid=j;
                                L=j;
                            end
                         end
                     end
                 end
                 S7(L).Cid='BS';
%                  if(sum(S7(S7(i).Cid).Cid)<n)
%                      LB7(S7(S7(i).Cid).Cid)=LB7(S7(S7(i).Cid).Cid)+1;
%                  end
             end
         end
     end
     end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%����ڵ����ϼ��ڵ�Cid�ľ���
DCid7 = zeros(1,B);
for i=1:1:n
    if (S7(i).E>0)
        if(sum(S7(i).Cid)~=149)
            DCid7(1,i)=sqrt((S7(i).xd-S7(S7(i).Cid).xd)^2+(S7(i).yd-S7(S7(i).Cid).yd)^2); 
        else
            DCid7(1,i)=S7(i).D;
        end
    end
end
for i=1:1:n
    if (DCid7(i)~=0 && DCid7(i)>=do)
        if(((PL+LB7(i))*rand(1,1)+(DCid7(i)/do)-1)>PL)
            Loss7(1,q)=Loss7(1,q)+1;
        end
    end
end 
 %%%%%%%%%%%%%%%%%%%%
 %�ܺ�ģ��
for i=1:1:n
    S6(i).EB=S6(i).E;
end
for i=1:1:n
 if(S7(i).E>0)
     if(sum(S7(i).Cid)~=149)
         d=sqrt((S7(i).xd-S7(S7(i).Cid).xd)^2 + (S7(i).yd-S7(S7(i).Cid).yd)^2);
        if (d>do)
                S7(i).E=S7(i).E-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB7(i) + Emp*4000*( d * d * d * d));
            end
            if (d<=do)
                S7(i).E=S7(i).E-((EPC+ETX)*(4000)+(ERX+EDA)*(4000)*LB7(i) + Efs*4000*( d * d)); 
            end
     else
          d=sqrt((S7(i).xd-BS.x)^2 + (S7(i).yd-BS.y)^2);
          if (d>do)
                S7(i).E=S7(i).E-(ETX*(4000)+(ERX+EDA)*(4000)*LB7(i) + Emp*4000*( d * d * d * d)); 
            end
            if (d<=do)
                S7(i).E=S7(i).E-(ETX*(4000)+(ERX+EDA)*(4000)*LB7(i) + Efs*4000*( d * d)); 
            end
     end
 end
end

if(countCHs7>0)%�д�ͷʱ����
    EB7=[];
    Eb7=1;
    Cb7=0;
    for i=1:1:n
        if (S7(i).type=='C')%���д�ͷ������������ʱ
            if (S7(i).EB>0 && S7(i).E>0)
                EB7(Eb7)=(S7(i).EB-S7(i).E)/S7(i).EB*100;
                if (EB7(Eb7)>Cb7)
                    Cb7=EB7(Eb7);
                end
                Eb7=Eb7+1;
            end
        end
    end
    if(length(EB7)>=1)%���д�ͷ��������������ʱ
        Ave7=sum(EB7)/length(EB7);
        STATISTICS.Ave7(r+1)=Ave7;
        STATISTICS.Cb7(r+1)=Cb7;
    else
        STATISTICS.Ave7(r+1)=0;
        STATISTICS.Cb7(r+1)=0;
    end
else
    STATISTICS.Ave7(r+1)=0;
    STATISTICS.Cb7(r+1)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ݴ����������ݽ�����Ϊ100%
for i=1:1:n
    if(S7(i).E>0)
       data7=data7+1; 
    end
end
end
STATISTICS.DATA3(q)=data3;
STATISTICS.DATA4(q)=data4;
STATISTICS.DATA5(q)=data5;
STATISTICS.DATA6(q)=data6;
STATISTICS.DATA7(q)=data7;
end
AvLoss3(1,o)=mean(Loss3);
AvLoss4(1,o)=mean(Loss4);
AvLoss5(1,o)=mean(Loss5);
AvLoss6(1,o)=mean(Loss6);
AvLoss7(1,o)=mean(Loss7);
VaLoss3(1,o)=std(Loss3);%���㶪ʧ���ݰ�
VaLoss4(1,o)=std(Loss4);
VaLoss5(1,o)=std(Loss5);
VaLoss6(1,o)=std(Loss6);
VaLoss7(1,o)=std(Loss7);
end
% q=1:10;
% figure(1)
% plot(q,STATISTICS.DATA3,'-b',q,STATISTICS.DATA4,'-r',q,STATISTICS.DATA6,'-g',q,STATISTICS.DATA5,'-y',q,STATISTICS.DATA7,'-k');
% legend('Balance','LEACH','LEACH-C','LEACH-N','LEACH-R');
% xlabel('X(Time)');
% ylabel('T(Data)');
% title('\bf Comparison of total data transfer per unit of energy for Balance, LEACH, LEACH-C, LEACH-N, LEACH-R');