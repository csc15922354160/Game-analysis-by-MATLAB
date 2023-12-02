# Game-analysis-by-MATLAB
use MATLAB to calculate and simulate the game, including the Stackelberg game and evolutionary game.
%% 斯塔克伯格博弈（内容包括：斯塔克伯格博弈的模型求解计算，海瑟矩阵，参数影响的二维和三维图）
%% 目标论文：《考虑保险理赔模式和政府补贴的再制造产品市场优化策略-斯塔克博格博弈模型》
% 模型R：求解结果
% 决策变量：再制造商：再制造零部件批发价格wRr;保险公司：保费费率折扣nR
% 利润函数，再制造商：PiRr=((wRr-cr+(1-s)*z)*ERr;保险公司：PiRb=(nR*P-wRr)*ERr+nR*P*(DRr-ERr);
% 保险产品规模：DRr=(d*(1-b)*(q-nR*P+s*z))/q ; 再制造零部件的期望需求：ERr=p*DRr
%&-------------------------------------------------------------------------结论1：
%第一步 根据上述定义参数
syms wRr nR cr s z ERr P DRr d b q p
%第二步 构建等式
DRr=(d*(1-b)*(q-nR*P+s*z))/q ;
ERr=p*DRr;
PiRr = (wRr-cr+(1-s)*z)*ERr;
PiRb = (nR*P-wRr)*ERr+nR*P*(DRr-ERr);
%第三步 求解最优决策（先保险公司Rb，再制造商Rr_optimized，再保险公司Rb_optimized,再产品保险需求DRr_optimized,再再制造零部件期望需求ERr_optimized, 再最大利润PiRr_optimized,PiRb_optimized）
sol_Rb = solve(diff(PiRb, nR) == 0 , nR);
disp('保险公司的反应函数为：'); pretty(simplify(sol_Rb));
PiRr_optimized = subs(PiRr, nR, sol_Rb);
sol_Rr = solve(diff(PiRr_optimized, wRr) == 0 , wRr);
disp('再制造商的最优批发价格为：'); pretty(simplify(sol_Rr));
sol_Rb_optimized = subs(sol_Rb, wRr, sol_Rr);
disp('保险商的最优费率折扣为：'); pretty(simplify(sol_Rb_optimized));
DRr_optimized = subs(DRr, [nR, wRr],[sol_Rb_optimized, sol_Rr]);
disp('产品保险的需求为：'); pretty(simplify(DRr_optimized ));
ERr_optimized = subs(ERr, DRr,DRr_optimized);
disp('再制造零部件的期望需求：'); pretty(simplify(ERr_optimized ));
PiRb_optimized = subs(PiRb, [nR, wRr],[sol_Rb_optimized, sol_Rr]);
disp('保险公司最优利润：'); pretty(simplify(PiRb_optimized));
PiRr_optimized1 = subs(PiRr_optimized,wRr,sol_Rr );
disp('再制造商最优利润：'); pretty(simplifyFraction(PiRr_optimized1));

%% ----------------------------------------------------------------------

%%画图，下面的为二维参数影响图
% R-1 p对wRr的影响  
syms wRr nR cr s z ERr P DRr d b q p
cr=0.4, P=1, d=1000, b=0.1, q=0.5;
s=1,z=0.1;
wRr=(q+p*cr+s*z-p*(1-s)*z)/(2*p);
wRr = ((2*p)/5 + 3/5)/(2*p);
s=0,z=0.1;
wRr=(q+p*cr+s*z-p*(1-s)*z)/(2*p)
wRr =((3*p)/10 + 1/2)/(2*p)
s=0,z=0;
wRr=(q+p*cr+s*z-p*(1-s)*z)/(2*p)
wRr =((2*p)/5 + 1/2)/(2*p)

figure('NumberTitle', 'off', 'Name', 's=1,z=1，s=0,z=1，s=0,z=0');
p=0.1:0.01:0.8;
%线1
wRr = ((2*p)/5 + 3/5)./(2*p)
plot(p,wRr, 'r.-')
hold on

%线2
wRr =((3*p)/10 + 1/2)./(2*p)
plot(p,wRr, 'r-')
hold on

%线3
wRr =((2*p)/5 + 1/2)./(2*p)
plot(p,wRr, 'b--')
legend('w^R_r','w^CR_r','w^RR_r')
hold on

xlabel('p');
ylabel('wr');
grid on
hold off

% R-2 p对再制造商利润的影响 
syms wRr nR cr s z ERr P DRr d b q p
cr=0.4, P=1, d=1000, b=0.1, q=0.5;
s=1,z=0.1;
PiRr=(d*(1-b)*(q-p*cr+s*z+p*(1-s)*z)^2)/(8*q)
PiRr =225*((2*p)/5 - 3/5)^2
s=0,z=0.1;
PiRr=(d*(1-b)*(q-p*cr+s*z+p*(1-s)*z)^2)/(8*q)
PiRr =225*((3*p)/10 - 1/2)^2
s=0,z=0;
PiRr=(d*(1-b)*(q-p*cr+s*z+p*(1-s)*z)^2)/(8*q)
PiRr =225*((2*p)/5 - 1/2)^2

figure('NumberTitle', 'off', 'Name', 's=1,z=1，s=0,z=1，s=0,z=0');
p=0.1:0.02:0.8;

PiRr =112.5.*((2.*p)/5 - 3/5).^2
plot(p,PiRr , 'r*-')

hold on

PiRr =112.5.*((3.*p)/10 - 1/2).^2
plot(p,PiRr , 'r-')

hold on

PiRr =112.5.*((2.*p)/5 - 1/2).^2
plot(p,PiRr , 'b^-')
legend('pi^R_r','pi^CR_r','pi^RR_r')
hold on

xlabel('p');
ylabel('pi');
grid on
hold off

%% %参数的影响图，三维的影响图
syms cn cr s z ERr P DRr d b q p u a 
cr=0.4, cn=0.8, z=0.1, d=1000, b=0.1, q=0.5;
s=1
u=s*z+p*(1-s)*z
a=2*(1-q)*(1-b*P);

w1= (2*a-(1-b)*(4*P*(1-q)-p*(2*q*cn+cr)+u))/(3*q*p*(1-b))
  w1 =(20*((27*p)/25 - 2*P + 191/100))/(27*p)
w2=(a-(1-b)*(2*P*(1-q)-p*(q*cn+2*cr)-s*z+2*p*(1-s)*z))/(3*p*(1-b))
  w2 =(10*((27*p)/25 - P + 109/100))/(27*p)

n=((1-b)*(2*P*(1+2*q)-p*(q*cn-cr)+5*s*z-p*(1-s)*z)-a)/(6*P*(1-b))
  n =(5*((37*P)/10 - 11/20))/(27*P)
 
figure('NumberTitle', 'off', 'Name', 's=1,z=1，s=0,z=1，s=0,z=0');
  [P,p]=meshgrid(0.2:0.01:1.1,0.2:0.01:0.8); 
  w1 =(20*((27*p)./25 - 2*P + 191/100))./(27*p)
  meshc(P,p,w1);title('meshc(P,p,w1)')   
hold on
xlabel('P');
ylabel('p');
zlabel('w1')
grid on
hold off

figure('NumberTitle', 'off', 'Name', 's=1,z=1,s=0,z=1，s=0,z=0');
  [P,p]=meshgrid(0.2:0.01:1.1,0.2:0.01:0.8); 
  w2 =(10*((27*p)./25 - P + 109/100))./(27*p)
  meshc(P,p,w2);title('meshc(P,p,w2)')   
hold on
xlabel('P');
ylabel('p');
zlabel('w2')
grid on
hold off

figure('NumberTitle', 'off', 'Name', 's=1,z=1，s=0,z=1，s=0,z=0');
  [cr,p]=meshgrid(0:0.01:1.1,0.2:0.01:0.8); 
  n =(5*p.*(cr - 2/5))./24 + 241/432
  meshc(cr,p,n);title('meshc(cr,p,n)')   
hold on
xlabel('cr');
ylabel('p');
zlabel('n')
grid on
hold off
