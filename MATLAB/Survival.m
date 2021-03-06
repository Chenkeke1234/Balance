clear
%1.初始参数设定模块
%.传感器节点区域界限(单位 M)
xm=100;
ym=100;
B=100;
rmax=4000;
%(1)汇聚节坐标给定
PI=[];
PI=0:pi/150:40*pi;
for i=1:1:4001
    Sink(i).x=xm*cos(PI(i))*rand(1,1);
    Sink(i).y=ym*sin(PI(i))*rand(1,1);
end
sink.x=0;
sink.y=0;
%区域内传器节数
n=100;
%簇头优化比例（当选簇头的概率）
P=0.05;
%能量模型（单位 焦）
%初始化能量模型
Eo=0.1;
%Eelec=Etx=Erx
ETX=50*0.000000001;%传输能量，每bit
ERX=50*0.000000001;%接收能量，每bit
%Transmit Amplifier types
Efs=10*0.000000000001;%自由空间损耗，每bit
Emp=0.0013*0.000000000001;%多径衰落损耗，每bit
%Data Aggregation Energy
EDA=5*0.000000001;%融合能耗，每bit
EPC=5*0.000000001;%每轮的感知收集单位信息的能耗，每bit
%最大循环次数

%算出参数 do，节点通信半径
do=sqrt(Efs/Emp)/2;
fo=ceil(do);
Et=0;
%数据传输变量
data3=0;
data4=0;
data5=0;
data6=0;
data7=0;
%2.无线传感器网络模型产生模块
%构建无线传感器网络,在区域内均匀投放100个节点,并画出图形
%初始化能级数组
e=zeros(4,B);%4表示区域数
a=zeros(4,B+1);
eav=zeros(4,B+1);
for i=1:1:n
    S3(i).xd=(rand(1,1)*2-1)*i;
    S3(i).yd=(rand(1,1)*2-1)*i;
%     S3(i).xd=(rand(1,1)*2-1)*i;
%     S3(i).yd=(rand(1,1)*2-1)*i;
    S4(i).xd=S3(i).xd;
    S4(i).yd=S3(i).yd;
    S5(i).xd=S3(i).xd;
    S5(i).yd=S3(i).yd;
    S6(i).xd=S3(i).xd;
    S6(i).yd=S3(i).yd;
    S7(i).xd=S3(i).xd;
    S7(i).yd=S3(i).yd;
    %计算每个节点与汇节点的距离
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
    %假设初始能量E0是经过一轮发送节点信息给汇节点后剩余的能量
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
    Et=Et+E3(i);%总能量
    %initially there are no cluster heads only nodes
    S3(i).type='N';%節點類型為普通
    S4(i).type='N';
    S5(i).type='N';
    S6(i).type='N';
    S7(i).type='N';
    S3(i).M=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%建立邻居节点表
L=10;%节点感知半径
R=50;%通信半径
for i=1:1:n
    Contigupus=zeros(1,n);%存储邻居节点编号
    Similar=zeros(1,n);
    A=zeros(1,n);
    for j=1:1:n
        d=sqrt((S3(i).xd-(S3(j).xd))^2 + (S3(i).yd-(S3(j).yd))^2);
        if(d<=2*L && d~=0 && S3(j).E>0)%节点间距离少于感知半径皆为邻居节点
            Contigupus(j)=j;   
        end
        if(d<L && d~=0 && S3(j).E>0)%挑选相似节点
            Similar(j)=j;   
        end
        if(S3(j).D>do && d<=do && d~=0 && S3(j).E>0 && S3(i).C==S3(j).C)
            A(j)=j;
        end
    end
    if (sum(find(Contigupus~=0))~=0)%判断矩阵是否为空
         S3(i).NB=find(Contigupus~=0);
         S7(i).NB=find(Contigupus~=0);
    else
         S3(i).NB='null';
         S7(i).NB='null';
    end
    if (sum(find(Similar~=0))~=0)%判断矩阵是否为空
         S3(i).S=find(Similar~=0);
    else
         S3(i).S='null';
    end
    if (sum(find(A~=0))~=0)%判断矩阵是否为空
         S3(i).a=length(find(A~=0))+1;
    else
         S3(i).a=1;
    end
end
d1=0.765*xm/2;
K=sqrt(0.5*n*do/pi)*xm/d1^2;
d2=xm/sqrt(2*pi*K);
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);
% S3(n+1).xd=sink.x;
% S3(n+1).yd=sink.y;
% S4(n+1).xd=sink.x;
% S4(n+1).yd=sink.y;
% S5(n+1).xd=sink.x;
% S5(n+1).yd=sink.y;
% S6(n+1).xd=sink.x;
% S6(n+1).yd=sink.y;
% S7(n+1).xd=sink.x;
% S7(n+1).yd=sink.y;
%死亡节点数
flag_first_dead3=0;
flag_teenth_dead3=0;
flag_all_dead3=0;
dead3=0;
first_dead3=0;
teenth_dead3=0;
all_dead3=0;
allive3=n;
%死亡节点数
allive4=n;
flag_first_dead4=0;
flag_teenth_dead4=0;
flag_all_dead4=0;
dead4=0;
first_dead4=0;
teenth_dead4=0;
all_dead4=0;
%死亡节点数
allive5=n;
flag_first_dead5=0;
flag_teenth_dead5=0;
flag_all_dead5=0;
dead5=0;
first_dead5=0;
teenth_dead5=0;
all_dead5=0;
%死亡节点数
allive6=n;
flag_first_dead6=0;
flag_teenth_dead6=0;
flag_all_dead6=0;
dead6=0;
first_dead6=0;
teenth_dead6=0;
all_dead6=0;
%死亡节点数
allive7=n;
flag_first_dead7=0;
flag_teenth_dead7=0;
flag_all_dead7=0;
dead7=0;
first_dead7=0;
teenth_dead7=0;
all_dead7=0;
 %%%%%%%%%%%%%%%%%%%%%
   %计算小层级能量
   for j=1:1:n
       e(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)+S3(j).E;%计算小层级总能量
       a(S3(j).C,S3(j).f)=a(S3(j).C,S3(j).f)+1;%计算该区域小层级的节点数
   end
   for j=1:1:n
       eav(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)./a(S3(j).C,S3(j).f);
       a(S3(j).C,B+1)= a(S3(j).C,B+1)+1;%计算各区域的节点个数
   end
Ea=Et/n;%Ea为剩余总平均能量
SCHCount=ceil(0.04*n);
N=zeros(SCHCount,n);%储存静态簇头筛选因子
R=zeros(1,n);%储存动态簇头筛选因子
JE=zeros(1,SCHCount);%用于判断区域是否存在动态簇头
countCHs3=0;%满足簇头从新筛选条件，簇头个数清零
for i=1:1:n
 N(S3(i).C,i)=(eav(S3(i).C,S3(i).f)*(S3(i).E-Ea))*(S3(i).a/n)/(4*Eo-S3(i).E)*(1-(S3(i).D/(1.42*B)));
%  if Ea>0
%      p(i)=P*S3(i).a/n;
%      if(S3(i).E>0)    
%            %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
%         N(S3(i).C,i)=(p(i)/(1-p(i)*mod(0,round(1/p(i))))*(S3(i).E/(3*eav(S3(i).C,S3(i).f))));%p(i)<1,大概率
%      end  
%   end 
end
[st,Sid]=max(N,[],2);%提取静态簇头筛选因子中的最大值及其ID
for j=1:1:SCHCount
    if (st(j)~=0 && a(S3(Sid(j)).C,B+1)>10) %满足条件则在区域建立静态簇头
       S3(Sid(j)).type='C';
       S3(Sid(j)).Cid='sink';%没有动态簇头时，静态簇头与汇节点sink连接
       countCHs3=countCHs3+1;
       if( S3(Sid(j)).f>fo )
           for z=1:1:n 
               if(S3(z).D<S3(Sid(j)).D && S3(z).C==S3(Sid(j)).C)%同一区域和到汇节点距离少于静态簇头时计算动态节点影响因子
                   R(z)=(4*(S3(Sid(j)).f-S3(z).f+1)*S3(z).f/((S3(Sid(j)).f+1)*(S3(Sid(j)).f+1)))*(1-S3(z).f/S3(Sid(j)).f)*(eav(S3(z).C,S3(z).f)*(S3(j).E-Ea))/(3*Eo-S3(z).E);
               end
           end
       [dy,Did]=max(R,[],2); %提取动态簇头筛选因子中的最大值及其ID  
       S3(Did).type='D';
       S3(Sid(j)).Cid=Did;%有动态簇头则重置为与动态簇头相连
       S3(Did).Cid='sink';%动态簇头与汇节点相连
       JE(S3(Did).C)=Did;
       countCHs3=countCHs3+1;
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%
%节点连接
LB=zeros(1,n);%接收信息计数器
for i=1:1:n
    if ( S3(i).type=='N' && S3(i).E>0 && S3(i).D<=do)
        S3(i).Cid='sink';
    elseif(S3(i).type=='N' && S3(i).E>0 && S3(i).D>do)
            if(JE(S3(i).C)~=0  && Sid(S3(i).C)~=0)%存在动、静态簇头
                d1=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                d2=sqrt((S3(i).xd-S3(JE(S3(i).C)).xd)^2+(S3(i).yd-S3(JE(S3(i).C)).yd)^2);
                if(d1<d2)
                    S3(i).Cid=Sid(S3(i).C);%保存节点所发送的特殊节点ID
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                    S3(i).Cid=JE(S3(i).C);
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                end
            elseif(Sid(S3(i).C)~=0 && JE(S3(i).C)==0)%不存在动态簇头，存在静态簇头
                d3=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                if(sum(S3(i).NB)==443 || d3<=do)%如果没有邻居节点，则直接与静态簇头相连,sum('null')==443
                     S3(i).Cid=Sid(S3(i).C);
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                elseif(sum(S3(i).NB)~=443 && d3>do)
                    b=[];
                    for z=1:1:length(S3(i).NB)
                        b(z)=S3(S3(i).NB(z)).E;
                    end
                    for j=1:1:length(S3(i).NB)
                        if(S3(S3(i).NB(j)).type~='N')%当邻居节点中有特殊节点时,连接特殊节点
                            S3(i).Cid=S3(i).NB(j);
                            LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            break;
                        else
                            if(S3(i).f>S3(S3(i).NB(j)).f && S3(S3(i).NB(j)).E==max(b))
                                 S3(i).Cid=S3(i).NB(j);%选择与层级少于自身并且能量最大的邻居节点相连
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                                 break;
                            else
                                 S3(i).Cid=Sid(S3(i).C);%否则直接与该域内静态簇头相连
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            end
                        end
                    end
                end
            else  %不存在簇头         
                if(S3(i).f>S3(S3(i).NB(j)).f)
                     S3(i).Cid=S3(i).NB(j);%选择与成绩少于自身的邻居节点相连
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                     S3(i).Cid='sink';%否则直接与该域内静态簇头相连
                end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%选择休眠等级
for i=1:1:n
    u=[];
    if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'或者为特殊节点
        for j=1:1:length(S3(i).S)
            u(j)=S3(S3(i).S(j)).M;
        end
        if(length(find(u==0))>0)
            if((S3(Sid(S3(i).C)).f>fo  && Sid(S3(i).C)~=0) || LB(i)==0)%存在动静态簇头的连接方式
                S3(i).M=2;
            else
                S3(i).M=1;
            end
        end
    end
end
%唤醒错误休眠节点（即一开始邻居节点是活跃的，后来休眠了）
for i=1:1:n
    u=[];
    if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'或者为特殊节点
        for j=1:1:length(S3(i).S)
            u(j)=S3(S3(i).S(j)).M;
        end
        if(length(find(u==0))==0)
            S3(i).M=0; 
        end
    end
end

for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r
    %死亡节点检查模块
dead3=0;
for i=1:1:n
    %检查有无死亡节点
    if (S3(i).E<=0)
        dead3=dead3+1; 
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead3==1)
           if(flag_first_dead3==0)
              first_dead3=r;
              flag_first_dead3=1;
           end
        end
        %10%的节点死亡时间
        if(dead3==0.05*n)
           if(flag_teenth_dead3==0)%每死亡5个节点
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
STATISTICS.DEAD3(r+1)=dead3;%死亡节点总数
STATISTICS.ALLIVE3(r+1)=allive3-dead3;%存活节点数
m=allive3-dead3;

%%%%%%%%%%%%%%%%%%%%%%%%%
%基站位置筛选
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
%更新邻居节点表
for i=1:1:n
    if(S3(i).E>0)
        Contigupus=zeros(1,n);%存储邻居节点编号
        Similar=zeros(1,n);
        A=zeros(1,n);
        for j=1:1:n
            d=sqrt((S3(i).xd-(S3(j).xd))^2 + (S3(i).yd-(S3(j).yd))^2);
            if(d<=2*L && d~=0 && S3(j).E>0)%节点间距离少于感知半径皆为邻居节点
                Contigupus(j)=j;   
            end
            if(d<L && d~=0 && S3(j).E>0)%挑选相似节点
                Similar(j)=j;   
            end
            if(S3(j).D>do && d<=do && d~=0 && S3(j).E>0 && S3(i).C==S3(j).C)
                A(j)=j;
            end
        end
        if (sum(find(Contigupus~=0))~=0)%判断矩阵是否为空
             S3(i).NB=find(Contigupus~=0);
        else
             S3(i).NB='null';
        end
        if (sum(find(Similar~=0))~=0)%判断矩阵是否为空
             S3(i).S=find(Similar~=0);
        else
             S3(i).S='null';
        end
        if (sum(find(A~=0))~=0)%判断矩阵是否为空
            S3(i).a=length(find(A~=0))+1;
        else
            S3(i).a=1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%当存活节点少于等于75时
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
   %计算小层级能量
e=zeros(4,B);%4表示区域数
a=zeros(4,B+1);
eav=zeros(4,B+1);
for j=1:1:n
   e(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)+S3(j).E;%计算小层级总能量
   a(S3(j).C,S3(j).f)=a(S3(j).C,S3(j).f)+1;%计算该区域小层级的节点数
end
for j=1:1:n
   eav(S3(j).C,S3(j).f)=e(S3(j).C,S3(j).f)./a(S3(j).C,S3(j).f);
   a(S3(j).C,B+1)= a(S3(j).C,B+1)+1;%计算各区域的节点个数
end
Et=0;
for i=1:1:n
    if(S3(i).E>0)
        Et=Et+S3(i).E;%总剩余能量
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ea=Et/(allive3-dead3);%Ea为剩余总平均能量
SCHCount=ceil(0.04*m);
N=zeros(SCHCount,n);%储存静态簇头筛选因子
R=zeros(1,n);%储存动态簇头筛选因子
JE=zeros(1,SCHCount);%用于判断区域是否存在动态簇头
countCHs3=0;%满足簇头从新筛选条件，簇头个数清零
if(m>=10)
for i=1:1:n
   N(S3(i).C,i)=(eav(S3(i).C,S3(i).f)*(S3(i).E-Ea))*(S3(i).a/n)/(4*Eo-S3(i).E)*(1-(S3(i).D/(1.42*B)));
%  if Ea>0
%      p(i)=P*S3(i).a/n;
%      if(S3(i).E>0)    
%            %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
%         N(S3(i).C,i)=(p(i)/(1-p(i)*mod(r,round(1/p(i))))*(S3(i).E/(3*eav(S3(i).C,S3(i).f))));%p(i)<1,大概率
%      end  
%   end 
end
[st,Sid]=max(N,[],2);%提取静态簇头筛选因子中的最大值及其ID
for j=1:1:SCHCount
if (st(j)~=0 && a(S3(Sid(j)).C,B+1)>10) %满足条件则在区域建立静态簇头
   S3(Sid(j)).type='C';
   S3(Sid(j)).Cid='sink';%没有动态簇头时，静态簇头与汇节点sink连接
   countCHs3=countCHs3+1;
   if( S3(Sid(j)).f>fo )
       for z=1:1:n 
           if(S3(z).D<S3(Sid(j)).D && S3(z).C==S3(Sid(j)).C)%同一区域和到汇节点距离少于静态簇头时计算动态节点影响因子
               R(z)=(4*(S3(Sid(j)).f-S3(z).f+1)*S3(z).f/((S3(Sid(j)).f+1)*(S3(Sid(j)).f+1)))*(1-S3(z).f/S3(Sid(j)).f)*(eav(S3(z).C,S3(z).f)*(S3(j).E-Ea))/(3*Eo-S3(z).E);
           end
       end
   [dy,Did]=max(R,[],2); %提取动态簇头筛选因子中的最大值及其ID  
   S3(Did).type='D';
   S3(Sid(j)).Cid=Did;%有动态簇头则重置为与动态簇头相连
   S3(Did).Cid='sink';%动态簇头与汇节点相连
   JE(S3(Did).C)=Did;
   countCHs3=countCHs3+1;
   end
end
end
end
STATISTICS.COUNTCHS3(r+1)=countCHs3;
%%%%%%%%%%%%%%%%%%%%%%
%节点连接
LB=zeros(1,n);%接收信息计数器
for i=1:1:n
    if ( S3(i).type=='N' && S3(i).E>0 && S3(i).D<=do)
        S3(i).Cid='sink';
    elseif(S3(i).type=='N' && S3(i).E>0 && S3(i).D>do)
            if(JE(S3(i).C)~=0  && Sid(S3(i).C)~=0)%存在动、静态簇头
                d1=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                d2=sqrt((S3(i).xd-S3(JE(S3(i).C)).xd)^2+(S3(i).yd-S3(JE(S3(i).C)).yd)^2);
                if(d1<d2)
                    S3(i).Cid=Sid(S3(i).C);%保存节点所发送的特殊节点ID
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                    S3(i).Cid=JE(S3(i).C);
                    LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                end
            elseif(Sid(S3(i).C)~=0 && JE(S3(i).C)==0)%不存在动态簇头，存在静态簇头
                d3=sqrt((S3(i).xd-S3(Sid(S3(i).C)).xd)^2+(S3(i).yd-S3(Sid(S3(i).C)).yd)^2);
                if(sum(S3(i).NB)==443 || d3<=do)%如果没有邻居节点，则直接与静态簇头相连,sum('null')==443
                     S3(i).Cid=Sid(S3(i).C);
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                elseif(sum(S3(i).NB)~=443 && d3>do)
                    b=[];
                    for z=1:1:length(S3(i).NB)
                        b(z)=S3(S3(i).NB(z)).E;
                    end
                    for j=1:1:length(S3(i).NB)
                        if(S3(S3(i).NB(j)).type~='N')%当邻居节点中有特殊节点时,连接特殊节点
                            S3(i).Cid=S3(i).NB(j);
                            LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            break;
                        else
                            if(S3(i).f>S3(S3(i).NB(j)).f && S3(S3(i).NB(j)).E==max(b))
                                 S3(i).Cid=S3(i).NB(j);%选择与层级少于自身并且能量最大的邻居节点相连
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                                 break;
                            else
                                 S3(i).Cid=Sid(S3(i).C);%否则直接与该域内静态簇头相连
                                 LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                            end
                        end
                    end
                end
            else  %不存在簇头         
                if(S3(i).f>S3(S3(i).NB(j)).f)
                     S3(i).Cid=S3(i).NB(j);%选择与成绩少于自身的邻居节点相连
                     LB(S3(i).Cid)=LB(S3(i).Cid)+1;
                else
                     S3(i).Cid='sink';%否则直接与该域内静态簇头相连
                end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
%选择休眠等级
for i=1:1:n
    if(S3(i).E>0)
        u=[];
        if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'或者为特殊节点
            for j=1:1:length(S3(i).S)
                u(j)=S3(S3(i).S(j)).M;
            end
            if(length(find(u==0))>0)
                if((S3(Sid(S3(i).C)).f>fo  && Sid(S3(i).C)~=0) || LB(i)==0)%存在动静态簇头的连接方式
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
%唤醒错误休眠节点（即一开始邻居节点是活跃的，后来休眠了）
for i=1:1:n
     if(S3(i).E>0)
        u=[];
        if(sum(S3(i).S)~=443 && S3(i).type=='N')%S3(i).S=='null'或者为特殊节点
            for j=1:1:length(S3(i).S)
                u(j)=S3(S3(i).S(j)).M;
            end
            if(length(find(u==0))==0)
                S3(i).M=0; 
            end
        end
        if(sum(S3(i).Cid)~=437 && S3(i).M~=0)%减去休眠节点的连接数
            LB(S3(i).Cid)=LB(S3(i).Cid)-1;
        end
     end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %为簇头补充能量
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %能耗模型
%节点每轮发送一次数据，每次为4000bit
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
if(countCHs3>0)%有簇头时计算
    EB3=[];
    Eb3=1;
    Cb3=0;
    for i=1:1:n
        if (S3(i).type=='C' || S3(i).type=='D')%当有簇头而不符合条件时
            if (S3(i).EB>0 && S3(i).E>0)
                EB3(Eb3)=(S3(i).EB-S3(i).E)/S3(i).EB*100;
                if (EB3(Eb3)>Cb3)
                    Cb3=EB3(Eb3);
                end
                Eb3=Eb3+1;
            end
        end
    end
    if(length(EB3)>=1)%当有簇头但都不符合条件时
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
% %重置休眠属
for i=1:1:n
    if(S3(i).E>0)
        S3(i).M=0;
    else
        S3(i).M=2;
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%leach算法
for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r
  %每过一个轮转周期(本程序为20次)使各节点的S（i）.G参数（该参数用于后面的簇选举，
  % 在该轮转周期内已当选过簇头的节点不能再当选）恢复为零，即50个周期，每周期20轮
  if(mod(r,round(1/P))==0)%round四舍五入取整函数, mod(x,y)为求余函数，y为除数，符号与除数相同
    for i=1:1:n
        S4(i).G=0;
    end
  end
  
  Et4=0;
  for i=1:1:n
      if(S4(i).E>0)
        Et4=S4(i).E+Et4;%El3为剩余总能量
      end
  end
%(2)死亡节点检查模块
dead4=0;
for i=1:1:n
    %检查有无死亡节点
    if (S4(i).E<=0)
        dead4=dead4+1; 
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead4==1)
           if(flag_first_dead4==0)
              first_dead4=r;
              flag_first_dead4=1;
           end
        end
        %10%的节点死亡时间
        if(dead4==0.1*n)
           if(flag_teenth_dead4==0)%如果满一个周期，20轮
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
STATISTICS.DEAD4(r+1)=dead4;%死亡节点总数
STATISTICS.ALLIVE4(r+1)=allive4-dead4;%存活节点数
m4=allive4-dead4;
Ea4=Et4/(allive4-dead4);%Et为总能量，Et*(1-r/rmax)为每一轮剩余的总能量，Ea为剩余总平均能量
%(4)簇头选举模块
countCHs4=0;
cluster4=1;
C4=[];
for i=1:1:n
%  if Ea4>0
%  p(i)=P*n*S4(i).E*S4(i).E/(Et*Ea4);% E3(i)= S3(i).E
 if(S4(i).E>0 && Ea4>0)
   temp_rand=rand;     
   if ( (S4(i).G)<=0)  
       %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
        if(temp_rand<= (P/(1-P*mod(r,round(1/P)))))%p(i)<1,大概率
            %p(i)/(1-p(i)*mod(r,round(1/p(i)))))为阈值计算公式
            countCHs4=countCHs4+1;
            S4(i).type='C';
            S4(i).Cid='sink';
            S4(i).G=round(1/P)-1;%p(i)越小，S3(i).G越大，大于0时不能竞选簇头
            C4(cluster4)=i;
            cluster4=cluster4+1;
        end      
    end
    % S3(i).G=S3(i).G-1;     
  end 
%  end
end
STATISTICS.COUNTCHS4(r+1)=countCHs4;%簇头个数

%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
LB4=zeros(1,n);
for i=1:1:n
   if ( S4(i).type=='N' && S4(i).E>0 )
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:cluster4-1
           temp=min(min_dis,sqrt((S4(i).xd-S4(C4(c)).xd)^2 + (S4(i).yd-S4(C4(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S4(i).Cid=C4(c);%保存普通所连接的簇头的编号
               LB4(C4(c))=LB4(C4(c))+1;%计算簇内节点个数
           end
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算数据传输量，数据接收率为100%
for i=1:1:n
    if(S4(i).E>0)
       data4=data4+1; 
    end
end
STATISTICS.DATA4(r+1)=data4;
%%%%%%%%%%%%%%%%%%%%%%%%
%能耗模型
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
if(countCHs4>0)%有簇头时计算
    EB4=[];
    Eb4=1;
    Cb4=0;
    for i=1:1:n
        if (S4(i).type=='C')%当有簇头而不符合条件时
            if (S4(i).EB>0 && S4(i).E>0)
                EB4(Eb4)=(S4(i).EB-S4(i).E)/S4(i).EB*100;
                if (EB4(Eb4)>Cb4)
                    Cb4=EB4(Eb4);
                end
                Eb4=Eb4+1;
            end
        end
    end
    if(length(EB4)>=1)%当有簇头但都不符合条件时
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
%LEACH-N算法
for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r
  %每过一个轮转周期(本程序为20次)使各节点的S（i）.G参数（该参数用于后面的簇选举，
  % 在该轮转周期内已当选过簇头的节点不能再当选）恢复为零，即50个周期，每周期20轮
  if(mod(r,round(1/P))==0)%round四舍五入取整函数, mod(x,y)为求余函数，y为除数，符号与除数相同
    for i=1:1:n
        S5(i).G=0;
    end
  end 
  Et5=0;
  for i=1:1:n
      if(S5(i).E>0)
        Et5=S5(i).E+Et5;%El3为剩余总能量
      end
  end
%(2)死亡节点检查模块
dead5=0;
for i=1:1:n
    %检查有无死亡节点
    if (S5(i).E<=0)
        dead5=dead5+1; 
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead5==1)
           if(flag_first_dead5==0)
              first_dead5=r;
              flag_first_dead5=1;
           end
        end
        %10%的节点死亡时间
        if(dead5==0.1*n)
           if(flag_teenth_dead5==0)%如果满一个周期，20轮
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
STATISTICS.DEAD5(r+1)=dead5;%死亡节点总数
STATISTICS.ALLIVE5(r+1)=allive5-dead5;%存活节点数
m5=allive5-dead5;
Ea5=Et5/(allive5-dead5);%Et为总能量，Et*(1-r/rmax)为每一轮剩余的总能量，Ea为剩余总平均能量

%%%%%%%%%%%%%%%%%%%%%%%%%
%更新邻居节点表
for i=1:1:n
    if(S5(i).E>0)
        A=zeros(1,n);%存储邻居节点编号
        for j=1:1:n
            d=sqrt((S5(i).xd-(S5(j).xd))^2 + (S5(i).yd-(S5(j).yd))^2);
            if(d<=do && d~=0 && S5(j).E>0)%节点间距离少于感知半径皆为邻居节点
                A(j)=j;   
            end
        end
        if (sum(find(A~=0))~=0)%判断矩阵是否为空
             S5(i).a=length(find(A~=0));
        else
             S5(i).a=0;
        end
    end
end

%(4)簇头选举模块
countCHs5=0;
cluster5=1;
C5=[];
for i=1:1:n
 if Ea5>0
 p(i)=P*S5(i).a/n;
 if(S5(i).E>0)
   temp_rand=rand;     
   if ( (S5(i).G)<=0)  
       %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))*(S5(i).E/(3*Eo))))%p(i)<1,大概率
            %p(i)/(1-p(i)*mod(r,round(1/p(i)))))为阈值计算公式
            countCHs5=countCHs5+1;
            S5(i).type='C';
            S5(i).Cid='sink';
            S5(i).G=round(1/p(i))-1;%p(i)越小，S3(i).G越大，大于0时不能竞选簇头
            C5(cluster5)=i;
            cluster5=cluster5+1;
        end      
    end
    % S3(i).G=S3(i).G-1;     
  end 
 end
end
STATISTICS.COUNTCHS5(r+1)=countCHs5;%簇头个数

%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
LB5=zeros(1,n);
for i=1:1:n
   if ( S5(i).type=='N' && S5(i).E>0 )
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:cluster5-1
           temp=min(min_dis,sqrt((S5(i).xd-S5(C5(c)).xd)^2 + (S5(i).yd-S5(C5(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S5(i).Cid=C5(c);%保存普通所连接的簇头的编号
               LB5(S5(i).Cid)=LB5(S5(i).Cid)+1;%计算簇内节点个数
           end
       end
%        if(min_dis<S5(i).D)
%            LB5(S5(i).Cid)=LB5(S5(i).Cid)+1;%计算簇内节点个数
%        else
%            S5(i).Cid='sink';
%        end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算数据传输量，数据接收率为100%
for i=1:1:n
    if(S5(i).E>0)
       data5=data5+1; 
    end
end
STATISTICS.DATA5(r+1)=data5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%能耗模型
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
if(countCHs5>0)%有簇头时计算
    EB5=[];
    Eb5=1;
    Cb5=0;
    for i=1:1:n
        if (S5(i).type=='C')%当有簇头而不符合条件时
            if (S5(i).EB>0 && S5(i).E>0)
                EB5(Eb5)=(S5(i).EB-S5(i).E)/S5(i).EB*100;
                if (EB5(Eb5)>Cb5)
                    Cb5=EB5(Eb5);
                end
                Eb5=Eb5+1;
            end
        end
    end
    if(length(EB5)>=1)%当有簇头但都不符合条件时
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
%LEACH-C算法
for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r
  %每过一个轮转周期(本程序为20次)使各节点的S（i）.G参数（该参数用于后面的簇选举，
  % 在该轮转周期内已当选过簇头的节点不能再当选）恢复为零，即50个周期，每周期20轮
  if(mod(r,round(1/P))==0)%round四舍五入取整函数, mod(x,y)为求余函数，y为除数，符号与除数相同
    for i=1:1:n
        S6(i).G=0;
    end
  end
  
  Et6=0;
for i=1:100
   if(S6(i).E>0)
        Et6=S6(i).E+Et6;%El3为剩余总能量
   end
end
%(2)死亡节点检查模块
dead6=0;
for i=1:1:n
    %检查有无死亡节点
    if (S6(i).E<=0)
        dead6=dead6+1; 
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead6==1)
           if(flag_first_dead6==0)
              first_dead6=r;
              flag_first_dead6=1;
           end
        end
        %10%的节点死亡时间
        if(dead6==0.1*n)
           if(flag_teenth_dead6==0)%如果满一个周期，20轮
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
STATISTICS.DEAD6(r+1)=dead6;%死亡节点总数
STATISTICS.ALLIVE6(r+1)=allive6-dead6;%存活节点数
m6=allive6-dead6;
Ea6=Et6/(allive6-dead6);%Et为总能量，Et*(1-r/rmax)为每一轮剩余的总能量，Ea为剩余总平均能量

%%%%%%%%%%%%%%%%%%%%%%%%
%当存活节点少于等于75时
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

%(4)簇头选举模块
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
STATISTICS.COUNTCHS6(r+1)=countCHs6;%簇头个数
%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
LB6=zeros(1,n);
for i=1:1:n
   if ( S6(i).type=='N' && S6(i).E>0 )
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:cluster6-1
           temp=min(min_dis,sqrt((S6(i).xd-S6(C6(c)).xd)^2 + (S6(i).yd-S6(C6(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S6(i).Cid=C6(c);%保存普通所连接的簇头的编号
               LB6(C6(c))=LB6(C6(c))+1;%计算簇内节点个数
           end
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算数据传输量，数据接收率为100%
for i=1:1:n
    if(S6(i).E>0)
       data6=data6+1; 
    end
end
STATISTICS.DATA6(r+1)=data6;

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
if(countCHs6>0)%有簇头时计算
    EB6=[];
    Eb6=1;
    Cb6=0;
    for i=1:1:n
        if (S6(i).type=='C')%当有簇头而不符合条件时
            if (S6(i).EB>0 && S6(i).E>0)
                EB6(Eb6)=(S6(i).EB-S6(i).E)/S6(i).EB*100;
                if (EB6(Eb6)>Cb6)
                    Cb6=EB6(Eb6);
                end
                Eb6=Eb6+1;
            end
        end
    end
    if(length(EB6)>=1)%当有簇头但都不符合条件时
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
%LEACH-R算法
for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环
    r 
BS.x=(rand(1,1)*2-1)*xm;
BS.y=(rand(1,1)*2-1)*ym;
for i=1:1:n
    S7(i).D=sqrt((S7(i).xd-BS.x)^2 + (S7(i).yd-BS.y)^2);
end
  %每过一个轮转周期(本程序为20次)使各节点的S（i）.G参数（该参数用于后面的簇选举，
  % 在该轮转周期内已当选过簇头的节点不能再当选）恢复为零，即50个周期，每周期20轮
  if(mod(r,round(1/P))==0)%round四舍五入取整函数, mod(x,y)为求余函数，y为除数，符号与除数相同
    for i=1:1:n
        S7(i).G=0;
    end
  end
  Et7=0;
  for i=1:1:n
      if(S7(i).E>0)
        Et7=S7(i).E+Et7;%El3为剩余总能量
      end
  end
%(2)死亡节点检查模块
dead7=0;
for i=1:1:n
    %检查有无死亡节点
    if (S7(i).E<=0)
        dead7=dead7+1; 
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead7==1)
           if(flag_first_dead7==0)
              first_dead7=r;
              flag_first_dead7=1;
           end
        end
        %10%的节点死亡时间
        if(dead7==0.1*n)
           if(flag_teenth_dead7==0)%如果满一个周期，20轮
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
STATISTICS.DEAD7(r+1)=dead7;%死亡节点总数
STATISTICS.ALLIVE7(r+1)=allive7-dead7;%存活节点数
m7=allive7-dead7;
Ea7=Et7/(allive7-dead7);%Et为总能量，Et*(1-r/rmax)为每一轮剩余的总能量，Ea为剩余总平均能量
%(4)簇头选举模块
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
           %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
            if(temp_rand<= (P/(1-P*mod(r,round(1/P)))))%p(i)<1,大概率
                %p(i)/(1-p(i)*mod(r,round(1/p(i)))))为阈值计算公式
                countCHs7=countCHs7+1;
                S7(i).type='C';
                S7(i).Cid='BS';
                S7(i).G=round(1/P)-1;%p(i)越小，S7(i).G越大，大于0时不能竞选簇头
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
           %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
            if (temp_rand <= (P/(1-P*mod(r,round(1/P))))*S7(i).Yn/(abs(Ymed)))%p(i)<1,大概率
                countCHs7=countCHs7+1;
                S7(i).type='C';
                S7(i).Cid='BS';
                S7(i).G=round(1/P)-1;%p(i)越小，S7(i).G越大，大于0时不能竞选簇头
                C7(cluster7)=i;
                cluster7=cluster7+1;
            end 
       else
           if (temp_rand <= (P/(1-P*mod(r,round(1/P))))*S7(i).Yn)%p(i)<1,大概率
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
STATISTICS.COUNTCHS7(r+1)=countCHs7;%簇头个数

%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
LB7=zeros(1,n);
ds=Inf;
Dm=0;

for i=1:1:n
   if(S7(i).E>0)
   if ( S7(i).type=='N')
       min_dis=Inf;%iInf表示正无穷大
       for c=1:1:cluster7-1
           temp=min(min_dis,sqrt((S7(i).xd-S7(C7(c)).xd)^2 + (S7(i).yd-S7(C7(c)).yd)^2));
           if ( temp<min_dis )
               min_dis=temp;
               S7(i).Cid=C7(c);%保存普通所连接的簇头的编号
               LB7(C7(c))=LB7(C7(c))+1;%计算簇内节点个数
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
             for j=1:1:n %计算与所有节点ds距离最小值
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
                 for j=1:1:n %找全场节点
                     if(J~=j)%S7(i).Cid~=j排除自己，还没排除自己的下级簇头
                         if(j~=i)%排除自己的下级簇头
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
 %%%%%%%%%%%%%%%%%%%%
 %能耗模型
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

if(countCHs7>0)%有簇头时计算
    EB7=[];
    Eb7=1;
    Cb7=0;
    for i=1:1:n
        if (S7(i).type=='C')%当有簇头而不符合条件时
            if (S7(i).EB>0 && S7(i).E>0)
                EB7(Eb7)=(S7(i).EB-S7(i).E)/S7(i).EB*100;
                if (EB7(Eb7)>Cb7)
                    Cb7=EB7(Eb7);
                end
                Eb7=Eb7+1;
            end
        end
    end
    if(length(EB7)>=1)%当有簇头但都不符合条件时
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
%计算数据传输量，数据接收率为100%
for i=1:1:n
    if(S7(i).E>0)
       data7=data7+1; 
    end
end
STATISTICS.DATA7(r+1)=data7;
end
r=0:4000;
STATISTICS.LIVE3=n-STATISTICS.DEAD3;
STATISTICS.LIVE4=n-STATISTICS.DEAD4;
STATISTICS.LIVE5=n-STATISTICS.DEAD5;
STATISTICS.LIVE6=n-STATISTICS.DEAD6;
STATISTICS.LIVE7=n-STATISTICS.DEAD7;
L3=[];
L4=[];
L5=[];
L6=[];
L7=[];
l=1;
for i=1:1:4000
    if(mod(i,round(500))==0)
       L3(l)=STATISTICS.LIVE3(i);
       L4(l)=STATISTICS.LIVE4(i);
       L5(l)=STATISTICS.LIVE5(i);
       L6(l)=STATISTICS.LIVE6(i);
       L7(l)=STATISTICS.LIVE7(i);
       l=l+1;
    end
end
h=1:1:8;
figure('name','节点存活情况比较')
plot(h,L3,'-b',h,L4,'-r',h,L6,'-g',h,L5,'-y',h,L7,'-k');
legend('Balance','LEACH','LEACH-C','LEACH-N','LEACH-R');
xlabel('x(time)');
ylabel('y(live)');
title('\bf Balance、LEACH、LEACH-C、LEACH-N、LEACH-R的存活节点数随时间的变化对比');

Z3=zeros(1,10);
ZM3=zeros(1,10);
Z4=zeros(1,10);
ZM4=zeros(1,10);
Z5=zeros(1,10);
ZM5=zeros(1,10);
Z6=zeros(1,10);
ZM6=zeros(1,10);
Z7=zeros(1,10);
ZM7=zeros(1,10);
for i=1:1:4001
    if(STATISTICS.Ave3(i)>=0 && STATISTICS.Ave3(i)<10)
        Z3(1)=Z3(1)+1;
        if(STATISTICS.Ave3(i)>ZM3(1))
            ZM3(1)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=10 && STATISTICS.Ave3(i)<20)
        Z3(2)=Z3(2)+1;
        if(STATISTICS.Ave3(i)>ZM3(2))
            ZM3(2)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=20 && STATISTICS.Ave3(i)<30)
        Z3(3)=Z3(3)+1;
        if(STATISTICS.Ave3(i)>ZM3(3))
            ZM3(3)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=30 && STATISTICS.Ave3(i)<40)
        Z3(4)=Z3(4)+1;
        if(STATISTICS.Ave3(i)>ZM3(4))
            ZM3(4)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=40 && STATISTICS.Ave3(i)<50)
        Z3(5)=Z3(5)+1;
        if(STATISTICS.Ave3(i)>ZM3(5))
            ZM3(5)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=50 && STATISTICS.Ave3(i)<60)
        Z3(6)=Z3(6)+1;
        if(STATISTICS.Ave3(i)>ZM3(6))
            ZM3(6)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=60 && STATISTICS.Ave3(i)<70)
        Z3(7)=Z3(7)+1;
        if(STATISTICS.Ave3(i)>ZM3(7))
            ZM3(7)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=70 && STATISTICS.Ave3(i)<80)
        Z3(8)=Z3(8)+1;
        if(STATISTICS.Ave3(i)>ZM3(8))
            ZM3(8)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=80 && STATISTICS.Ave3(i)<90)
        Z3(9)=Z3(9)+1;
        if(STATISTICS.Ave3(i)>ZM3(9))
            ZM3(9)=STATISTICS.Ave3(i);
        end
    elseif(STATISTICS.Ave3(i)>=90 && STATISTICS.Ave3(i)<=100)
        Z3(10)=Z3(10)+1;
        if(STATISTICS.Ave3(i)>ZM3(10))
            ZM3(10)=STATISTICS.Ave3(i);
        end
    end
end

for i=1:1:4001
    if(STATISTICS.Ave4(i)>=0 && STATISTICS.Ave4(i)<10)
        Z4(1)=Z4(1)+1;
        if(STATISTICS.Ave4(i)>ZM4(1))
            ZM4(1)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=10 && STATISTICS.Ave4(i)<20)
        Z4(2)=Z4(2)+1;
        if(STATISTICS.Ave4(i)>ZM4(2))
            ZM4(2)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=20 && STATISTICS.Ave4(i)<30)
        Z4(3)=Z4(3)+1;
        if(STATISTICS.Ave4(i)>ZM4(3))
            ZM4(3)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=30 && STATISTICS.Ave4(i)<40)
        Z4(4)=Z4(4)+1;
        if(STATISTICS.Ave4(i)>ZM4(4))
            ZM4(4)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=40 && STATISTICS.Ave4(i)<50)
        Z4(5)=Z4(5)+1;
        if(STATISTICS.Ave4(i)>ZM4(5))
            ZM4(5)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=50 && STATISTICS.Ave4(i)<60)
        Z4(6)=Z4(6)+1;
        if(STATISTICS.Ave4(i)>ZM4(6))
            ZM4(6)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=60 && STATISTICS.Ave4(i)<70)
        Z4(7)=Z4(7)+1;
        if(STATISTICS.Ave4(i)>ZM4(7))
            ZM4(7)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=70 && STATISTICS.Ave4(i)<80)
        Z4(8)=Z4(8)+1;
        if(STATISTICS.Ave4(i)>ZM4(8))
            ZM4(8)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=80 && STATISTICS.Ave4(i)<90)
        Z4(9)=Z4(9)+1;
        if(STATISTICS.Ave4(i)>ZM4(9))
            ZM4(9)=STATISTICS.Ave4(i);
        end
    elseif(STATISTICS.Ave4(i)>=90 && STATISTICS.Ave4(i)<=100)
        Z4(10)=Z4(10)+1;
        if(STATISTICS.Ave4(i)>ZM4(10))
            ZM4(10)=STATISTICS.Ave4(i);
        end
    end
end

for i=1:1:4001
    if(STATISTICS.Ave5(i)>=0 && STATISTICS.Ave5(i)<10)
        Z5(1)=Z5(1)+1;
        if(STATISTICS.Ave5(i)>ZM5(1))
            ZM5(1)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=10 && STATISTICS.Ave5(i)<20)
        Z5(2)=Z5(2)+1;
        if(STATISTICS.Ave5(i)>ZM5(2))
            ZM5(2)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=20 && STATISTICS.Ave5(i)<30)
        Z5(3)=Z5(3)+1;
        if(STATISTICS.Ave5(i)>ZM5(3))
            ZM5(3)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=30 && STATISTICS.Ave5(i)<40)
        Z5(4)=Z5(4)+1;
        if(STATISTICS.Ave5(i)>ZM5(4))
            ZM5(4)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=40 && STATISTICS.Ave5(i)<50)
        Z5(5)=Z5(5)+1;
        if(STATISTICS.Ave5(i)>ZM5(5))
            ZM5(5)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=50 && STATISTICS.Ave5(i)<60)
        Z5(6)=Z5(6)+1;
        if(STATISTICS.Ave5(i)>ZM5(6))
            ZM5(6)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=60 && STATISTICS.Ave5(i)<70)
        Z5(7)=Z5(7)+1;
        if(STATISTICS.Ave5(i)>ZM5(7))
            ZM5(7)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=70 && STATISTICS.Ave5(i)<80)
        Z5(8)=Z5(8)+1;
        if(STATISTICS.Ave5(i)>ZM5(8))
            ZM5(8)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=80 && STATISTICS.Ave5(i)<90)
        Z5(9)=Z5(9)+1;
        if(STATISTICS.Ave5(i)>ZM5(9))
            ZM5(9)=STATISTICS.Ave5(i);
        end
    elseif(STATISTICS.Ave5(i)>=90 && STATISTICS.Ave5(i)<=100)
        Z5(10)=Z5(10)+1;
        if(STATISTICS.Ave5(i)>ZM5(10))
            ZM5(10)=STATISTICS.Ave5(i);
        end
    end
end

for i=1:1:4001
    if(STATISTICS.Ave6(i)>=0 && STATISTICS.Ave6(i)<10)
        Z6(1)=Z6(1)+1;
        if(STATISTICS.Ave6(i)>ZM6(1))
            ZM6(1)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=10 && STATISTICS.Ave6(i)<20)
        Z6(2)=Z6(2)+1;
        if(STATISTICS.Ave6(i)>ZM6(2))
            ZM6(2)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=20 && STATISTICS.Ave6(i)<30)
        Z6(3)=Z6(3)+1;
        if(STATISTICS.Ave6(i)>ZM6(3))
            ZM6(3)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=30 && STATISTICS.Ave6(i)<40)
        Z6(4)=Z6(4)+1;
        if(STATISTICS.Ave6(i)>ZM6(4))
            ZM6(4)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=40 && STATISTICS.Ave6(i)<50)
        Z6(5)=Z6(5)+1;
        if(STATISTICS.Ave6(i)>ZM6(5))
            ZM6(5)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=50 && STATISTICS.Ave6(i)<60)
        Z6(6)=Z6(6)+1;
        if(STATISTICS.Ave6(i)>ZM6(6))
            ZM6(6)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=60 && STATISTICS.Ave6(i)<70)
        Z6(7)=Z6(7)+1;
        if(STATISTICS.Ave6(i)>ZM6(7))
            ZM6(7)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=70 && STATISTICS.Ave6(i)<80)
        Z6(8)=Z6(8)+1;
        if(STATISTICS.Ave6(i)>ZM6(8))
            ZM6(8)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=80 && STATISTICS.Ave6(i)<90)
        Z6(9)=Z6(9)+1;
        if(STATISTICS.Ave6(i)>ZM6(9))
            ZM6(9)=STATISTICS.Ave6(i);
        end
    elseif(STATISTICS.Ave6(i)>=90 && STATISTICS.Ave6(i)<=100)
        Z6(10)=Z6(10)+1;
        if(STATISTICS.Ave6(i)>ZM6(10))
            ZM6(10)=STATISTICS.Ave6(i);
        end
    end
end

for i=1:1:4001
    if(STATISTICS.Ave7(i)>=0 && STATISTICS.Ave7(i)<10)
        Z7(1)=Z7(1)+1;
        if(STATISTICS.Ave7(i)>ZM7(1))
            ZM7(1)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=10 && STATISTICS.Ave7(i)<20)
        Z7(2)=Z7(2)+1;
        if(STATISTICS.Ave7(i)>ZM7(2))
            ZM7(2)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=20 && STATISTICS.Ave7(i)<30)
        Z7(3)=Z7(3)+1;
        if(STATISTICS.Ave7(i)>ZM7(3))
            ZM7(3)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=30 && STATISTICS.Ave7(i)<40)
        Z7(4)=Z7(4)+1;
        if(STATISTICS.Ave7(i)>ZM7(4))
            ZM7(4)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=40 && STATISTICS.Ave7(i)<50)
        Z7(5)=Z7(5)+1;
        if(STATISTICS.Ave7(i)>ZM7(5))
            ZM7(5)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=50 && STATISTICS.Ave7(i)<60)
        Z7(6)=Z7(6)+1;
        if(STATISTICS.Ave7(i)>ZM7(6))
            ZM7(6)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=60 && STATISTICS.Ave7(i)<70)
        Z7(7)=Z7(7)+1;
        if(STATISTICS.Ave7(i)>ZM7(7))
            ZM7(7)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=70 && STATISTICS.Ave7(i)<80)
        Z7(8)=Z7(8)+1;
        if(STATISTICS.Ave7(i)>ZM7(8))
            ZM7(8)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=80 && STATISTICS.Ave7(i)<90)
        Z7(9)=Z7(9)+1;
        if(STATISTICS.Ave7(i)>ZM7(9))
            ZM7(9)=STATISTICS.Ave7(i);
        end
    elseif(STATISTICS.Ave7(i)>=90 && STATISTICS.Ave7(i)<=100)
        Z7(10)=Z7(10)+1;
        if(STATISTICS.Ave7(i)>ZM7(10))
            ZM7(10)=STATISTICS.Ave7(i);
        end
    end
end

z=1:1:10;
figure('name','簇头能耗分析')
plot(z,Z3,'-b',z,Z4,'-r',z,Z5,'-g',z,Z6,'-y',z,Z7,'-k');
legend('Balance','LEACH','LEACH-C','LEACH-N','LEACH-R');
xlabel('x(能耗区间)');
ylabel('y(簇头统计总数)');
title('\bf Balance、LEACH、LEACH-C、LEACH-N、LEACH-R的簇头能耗分析比较');
hold on;
figure('name','簇头能耗分析')
plot(z,ZM3,'-b',z,ZM4,'-r',z,ZM5,'-g',z,ZM6,'-y',z,ZM7,'-k');
legend('Balance','LEACH','LEACH-C','LEACH-N','LEACH-R');
xlabel('x(能耗区间)');
ylabel('y(能耗率)');
title('\bf Balance、LEACH、LEACH-C、LEACH-N、LEACH-R的各区间簇头能耗最大值分布');
y3=[];
y4=[];
y5=[];
y6=[];
y7=[];
t3=1;
t4=1;
t5=1;
t6=1;
t7=1;
for i=1:1:4001
    if(mod(i,round(1/0.00125))==0 && i<=4000)
        y3(t3)=STATISTICS.COUNTCHS3(i);
        t3=t3+1;
    end
    if(mod(i,round(1/0.00125))==0 && i<=4000)
        y4(t4)=STATISTICS.COUNTCHS4(i);
        t4=t4+1;
    end
    if(mod(i,round(1/0.00125))==0 && i<=4000)
        y5(t5)=STATISTICS.COUNTCHS5(i);
        t5=t5+1;
    end
    if(mod(i,round(1/0.00125))==0 && i<=4000)
        y6(t6)=STATISTICS.COUNTCHS6(i);
        t6=t6+1;
    end
    if(mod(i,round(1/0.00125))==0 && i<=4000)
        y7(t7)=STATISTICS.COUNTCHS7(i);
        t7=t7+1;
    end
end
figure('name','Balance簇头数目')
x3=[800 1600 2400 3200 4000];
bar(x3,y3,0.5);
legend('Balance');
xlabel('x(轮数)');
ylabel('y(Count)');
% set(gca,'YTickLabel',{'0','2','4','6','8','10'}) 
title('\bf Balance簇头数随时间的变化');

figure('name','LEACH簇头数目')
x4=[800 1600 2400 3200 4000];
bar(x4,y4,0.5);
legend('LEACH');
xlabel('x(轮数)');
ylabel('y(Count)');
title('\bf LEACH簇头数随时间的变化');

figure('name','LEACH-N簇头数目')
x5=[800 1600 2400 3200 4000];
bar(x5,y5,0.5);
legend('LEACH-N');
xlabel('x(轮数)');
ylabel('y(Count)');
title('\bf LEACH-N簇头数随时间的变化');

figure('name','LEACH-C簇头数目')
x6=[800 1600 2400 3200 4000];
bar(x6,y6,0.5);
legend('LEACH-C');
xlabel('x(轮数)');
ylabel('y(Count)');
title('\bf LEACH-C簇头数随时间的变化');

figure('name','LEACH-R簇头数目')
x7=[800 1600 2400 3200 4000];
bar(x7,y7,0.5);
legend('LEACH-R');
xlabel('x(轮数)');
ylabel('y(Count)');
title('\bf LEACH-R簇头数随时间的变化');
