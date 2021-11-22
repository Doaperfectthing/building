 

%% ==================================������ʼ��=====================================
clc,clear;
%ʱ������
N =24;
%ʱ�䲽��
dt=1;
%ʱ������
t=1:dt:24;
%% �������
%̫�������
ss=[161,79,290,112,494,538,608,504,528,566,495,527,549,485,363,322,243,259,279,215,91,63,39,30];
%�ļ������¶ȣ����϶ȣ�
th0=[22.20,22.80,23.70,23.00,22.10,24.45,33.31,31.64,28.44,39.51,37.18,40.97,42.45,32.67,25.08,24.47,25.20,25.20,25.10,25.30,25.40,25.00,24.80,24.50];
%���������¶ȣ����϶ȣ�
tc0=[2.20,2.80,3.70,3.00,2.10,4.45,13.31,11.64,8.44,9.51,7.18,10.97,12.45,12.67,5.08,4.47,5.20,5.20,5.10,5.30,5.40,5.00,4.80,4.50];
%�����ڲ������䶯
Qint=[150.88,119.78,118.35,94.57,82.84,95.48,87.02,103.91,144.29,154.83,300.01,276.20,357.53,294.79,354.58,365.04,385.57,341.54,325.83,278.13,203.68,198.42,182.04,161.52];
%���ڷ����¶ȣ����϶ȣ�
td=[22.20,22.80,23.70,23.00,22.10,23.60,24.00,23.20,24.30,25.80,26.30,27.50,27.70,26.70,25.70,25.70,25.20,25.20,25.10,25.30,25.40,25.00,24.80,24.50];
%�¶����ʶȷ�Χ
tb=[24  24	24	24	24	23	22	21	21	21	21	21	21	21	21	21	21	21  22	23	24	24	24	24;
    17	17	17	17	17	18	19	20	20	20	20	20	20	20	20	20	20	20	19	18	17	17	17	17];

%% �¶ȵ�λת��
%�����¶ȣ�ת���϶ȣ�
tt=tc0;
ttk=tt+273.15;
t0=ttk;        
%���ڷ����¶ȣ�ת���϶ȣ�
tdk=td+273.15;
%�¶����ʶȷ�Χ��ת���϶ȣ�
tbk=tb+273.15;
%�����ڲ������䶯
qint=Qint;
%% ����¥�����
 
a12=27.5;%ǽ�����m2
a13=23;  %ǽ�����m2
a14=27.5;%ǽ�����m2
a15=20.2;%ǽ�����m2
aw15=2.8;%�������m2
m=1000;
s=3600;
c12=7.9e5/m;  %ǽ������kj/k  
c13=6.6e5/m;  %ǽ������kj/k
c14=7.9e5/m;  %ǽ������kj/k
c15=2.6e7/m;  %ǽ������kj/k
cr=2.8e5/m;   %���ڿ�������kj/k
r12=0.0640*m;   %ǽ������k/kw  
r13=0.0768*m;   %ǽ������k/kw
r14=0.0640*m;   %ǽ������k/kw
r15=0.0299*m;   %ǽ������k/kw
rw15=0.0868*m;  %��������k/kw

h0=298.6e3/m;     %���� kj/kg
s0=6.86e3/m;      %���� kj/kgk
cp=1005/m;   %��ѹ������kj/kgk  
cv=718/m;    %�����������kj/kgk
gr=287/m;    %���峣��kj/kgk

tw=0.9;   %����������
xs=0.4;   %����������
mr=8*4*3*1.205; %���ڿ����������� 

%% ���ʶȲ���
pen=50;         %���ʶȳͷ�1 kwh/k
pex=3.8e8;      %���ʶȳͷ�2 kwh/k
ft=1;           %���ʶ��ɳ�ϵ��
  
%% �ȱò���
ls=50;
md=0.52*3600/ls;    %������������kg/h  *3600/300
mn=0.30*3600/ls;    %������������kg/h   *3600/300
COP=3.2;
%�������
k=711/(1000*3600^3);%kW h3/kg3   *3600^3
 
ps=135; %pa
ad=1.29;%kg/m3
eff=0.15;

%�ͷ��¶�
dpmax=3;          %�����¶�k
pmax=30+273.15;   %�¶�����k
pmin=15+273.15;   %�¶�����k

%�������� 
% dmmax=0.2*md;%��������
% mmax=md;     %��������
% mmin=0;      %��������
% 
% mMax =  [mmax; 0; 0 ;0 ;0];
% mMin =  [mmin; 0; 0 ;0 ;0];
% dmMax = [dmmax];



 
%% ====================================�Ż�ģ��===================================== 
%
% Import the problem dimensions.
 
nx = 5;
nu = 1;
nw = 4;
  
piOLOC=sdpvar(nu,N);
xOLOC =sdpvar(nx,N+1);
%   miOLOC=sdpvar(nx,N);  
  
%% ����ģ��Լ������
miOLOC = zeros(1,N);
for i=1:24
    if 5<=i && i<=18
miOLOC(:,i) = md;
    else
miOLOC(:,i) = mn;
    end
end 



   
% ������̬����.
% A = zeros(nx,nx,N);
% ad=[0.8715 0.0033 0.0028 0.0033 0.0062;
%     0.0012 0.9976 0      0      0;
%     0.0012 0      0.9976 0      0;
%     0.0012 0      0      0.9976 0;
%     0.0001 0      0      0      0.9998];
% an=[0.9183 0.0033 0.0028 0.0033 0.0062;
%     0.0012 0.9976 0      0      0;
%     0.0012 0      0.9976 0      0;
%     0.0012 0      0      0.9976 0;
%     0.0001 0      0      0      0.9998];
% for i=1:24
%     if 5<=i && i<=18
% A(:,:,i) = ad;
%     else
% A(:,:,i) = an;
%     end
% end 
  
% ���� nx x nu ���ƾ���, B(k).
% B = zeros(nx,nu,N);
% bd=[0.1106 0 0 0 0]';
% bn=[0.0638 0 0 0 0]';

% for i=1:24
%     if 5<=i && i<=18
% B(:,:,i) = bd;
%     else
% B(:,:,i) = bn;
%     end
% end 


% ���� nx x nw ���ž���, G.
% E = zeros(nx,nw);
% E = [ 0      0      0      0.0023;
%       0.0012 0      0      0;
%       0      0.0012 0      0;
%       0      0      0.0012 0;
%       0      0      0      0.0001];
% w=zeros(nw,N);
% w=zeros(5,24);
% w(1,:) =(1/cr)*(tw*aw15*ss/1000+qint);
% w(2,:) =tdk/(c12*r12);
% w(3,:) =tdk/(c13*r13);
% w(4,:) =tdk/(c14*r14);
% w(5,:) =tdk/(c15*r15)+(1/c15)*xs*aw15*ss/1000;
  
% xOLOC(:,1) = x0;
% for kp = 1 : N
% xOLOC(:,kp+1) = A(:,:,kp)*xOLOC(:,kp) + B(:,:,kp)*piOLOC(:,kp) + w(:,kp);
%     
% end
  

%     F= F+(repmat(mMin,1,N) <= miOLOC <= repmat(mMax,1,N));
%     F=F+  (abs(miOLOC(1,kp) - miOLOC(1,kp+1))<= dmMax);



   
   % ״̬Լ��
    x0=18*ones(nx,1)+273.15; %�����ʼ�¶�
    qint=-qint/2; %�为�� 
%     qint=-qint;   %�为��
    xOLOC(:,1) = x0;
    
    for kp = 1 : N
    a=zeros(5,5);
    a(1,1)=(-1/cr)*(1/r12+1/r13+1/r14+1/r15+1/rw15);
    a(1,2)=(1/cr)*(1/r12);
    a(1,3)=(1/cr)*(1/r13);
    a(1,4)=(1/cr)*(1/r14);
    a(1,5)=(1/cr)*(1/r15+1/rw15);
    a(2,1)=(1/c12)*(1/r12);
    a(3,1)=(1/c13)*(1/r13);
    a(4,1)=(1/c14)*(1/r14);
    a(5,1)=(1/c15)*(1/r15);
    a(2,2)=(-2/c12)*(1/r12);
    a(3,3)=(-2/c13)*(1/r13);
    a(4,4)=(-2/c14)*(1/r14);
    a(5,5)=(-2/c15)*(1/r15);
    a=a+eye(5,5);
    
    b=(1/cr)*cp*miOLOC(1,kp)*(piOLOC(1,kp)-xOLOC(1,kp));
    
    w=zeros(5,24);
    w(1,:) =(1/cr)*(tw*aw15*ss/1000+qint);
    w(2,:) =tdk/(c12*r12);
    w(3,:) =tdk/(c13*r13);
    w(4,:) =tdk/(c14*r14);
    w(5,:) =tdk/(c15*r15)+(1/c15)*xs*aw15*ss/1000;
        
    xOLOC(:,kp+1) = a*xOLOC(:,kp) + b*[1;0;0;0;0] + w(:,kp);
    end
% ���Ʊ���Լ��
F= (pmin <= piOLOC <= pmax);
for kp = 1 : N-1
F=F+  (abs(piOLOC(1,kp) - piOLOC(1,kp+1))<= dpmax);
end
% ���ʶ�Լ��
F=F+  (tbk(2,:)<= xOLOC(1,2:N+1)<= tbk(1,:));
F=F+ (piOLOC>=xOLOC(1,2:N+1));      
 
%%                                              Define the cost function.
%�������
% p1=sum((1-t0(1,:)./xOLOC(1,2:N+1)).*((tdk(1,:)-xOLOC(1,2:N+1))/r12) +(tdk(1,:)-xOLOC(1,2:N+1))/r13+(tdk(1,:)-xOLOC(1,2:N+1))/r14+(tdk(1,:)-xOLOC(1,2:N+1))/r15);
% p2=sum(miOLOC(1,:).*(cp*abs(piOLOC(1,:)-xOLOC(1,1:N))-cv*t0(1,:).*log(piOLOC(1,:)./xOLOC(1,1:N))));
% p1=sum(t0(1,:).*(  (   (tdk(1,:)-xOLOC(1,2:N+1)).^2   )./(   r12*tdk(1,:).*xOLOC(1,2:N+1)  ) + ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r13*tdk(1,:).*xOLOC(1,2:N+1)) + ...
%                  ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r14*tdk(1,:).*xOLOC(1,2:N+1)) + ((t0(1,:)-xOLOC(1,2:N+1)).^2)./(r15*t0(1,:).*xOLOC(1,2:N+1))));
% p2=sum(t0(1,:).*miOLOC(1,:)*cp.*( (piOLOC(1,:)-xOLOC(1,1:N)).^2 ) ./ (piOLOC(1,:).*xOLOC(1,1:N)) );
% p1=sum(t0(1,:).*(  (   (tdk(1,:)-xOLOC(1,2:N+1))/r12 ).*(   1./tdk(1,:) - 1./xOLOC(1,2:N+1)  ) + ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r13*tdk(1,:).*xOLOC(1,2:N+1)) + ...
%                  ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r14*tdk(1,:).*xOLOC(1,2:N+1)) + ((t0(1,:)-xOLOC(1,2:N+1)).^2)./(r15*t0(1,:).*xOLOC(1,2:N+1))));
% p2=sum(t0(1,:).*miOLOC(1,:)*cp.* (piOLOC(1,:)-xOLOC(1,1:N)) .*(piOLOC(1,:)-xOLOC(1,1:N))./ (piOLOC(1,:).*xOLOC(1,1:N))) ;
%������������Ļ������
% p1=sum(((1-t0./xOLOC(1,1:N)).*(pos(xOLOC(1,1:N)-tdk)/r12) + pos(xOLOC(1,1:N)-tdk)/r13 + pos(xOLOC(1,1:N)-tdk)/r14 + (xOLOC(1,1:N)-t0)/r15));
p1=sum(t0(1,:).*(  (   (tdk(1,:)-xOLOC(1,2:N+1)).^2   )./(   r12*tdk(1,:).*xOLOC(1,2:N+1)  ) + ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r13*tdk(1,:).*xOLOC(1,2:N+1)) + ...
                 ((tdk(1,:)-xOLOC(1,2:N+1)).^2)./(r14*tdk(1,:).*xOLOC(1,2:N+1)) + ((t0(1,:)-xOLOC(1,2:N+1)).^2)./(r15*t0(1,:).*xOLOC(1,2:N+1))));

%���������������Ļ������
ppp=log(piOLOC)-log(xOLOC(1,1:N));
p2=sum(miOLOC.*(cp*(piOLOC-xOLOC(1,1:N))-cv*t0.*ppp));
%����״̬�仯����Ļ������
pp=log(xOLOC(1,1:N))-log(xOLOC(1,2:N+1));
p3=sum(mr*(cp*pos(xOLOC(1,1:N)-xOLOC(1,2:N+1))-cv*t0.*(pp)));
% p2=norm((cp.*miOLOC(1,:).*(piOLOC(1,:)-xOLOC(1,1:N)-t0(1,:).*(piOLOC(1,:)./xOLOC(1,1:N)-1))),1);
% p8=norm(((1-t0(1,:)./xOLOC(1,1:N)).*cp.*mr.*(xOLOC(1,2:N+1)-xOLOC(1,1:N))),1);

% delapx(1,:)=((1-t0(1,:)./xOLOC(1,1:N)).*(piOLOC(1,:)-xOLOC(1,1:N)));
% p2=sum((cp*miOLOC(1,1).*delapx));
% p3=sum(mr*(cp*(xOLOC(1,2:N+1)-xOLOC(1,1:N))-cv*t0(1,1:N).*log(xOLOC(1,2:N+1)./xOLOC(1,1:N))));
% delapx(1,:)=(piOLOC(1,:)-xOLOC(1,1:N)).*((1+t0(1,:)./xOLOC(1,1:N))); 
% delapx(1,:)=(piOLOC(1,:)-xOLOC(1,1:N))-t0(1,:).*log(piOLOC(1,:)./xOLOC(1,1:N));  %px
% p7=norm((cp*miOLOC(1,1).*delapx),1);
%     
%   pm=sum(mr*cp*(xOLOC(1,2:N+1)-xOLOC(1,1:N)));
%   pn=sum(mr*t0(1,:)*cv*0.00001);
%   md=md/1200;
%   pc=sum((cp*miOLOC(1,:).*((1-t0(1,:)./piOLOC(1,:))-((1-t0(1,:)./xOLOC(1,1:N))))));
px=p1+p2-p3;
% px=p1+p2+p3;

%�����
% p4=sum(((tdk(1,:)-xOLOC(1,1:N))/r12) +(tdk(1,:)-xOLOC(1,1:N))/r13+(tdk(1,:)-xOLOC(1,1:N))/r14+(t0(1,:)-xOLOC(1,1:N))/r15);
p5=sum(miOLOC.*(cp.*(piOLOC-xOLOC(1,1:N))));
% p9=sum(cp.*mr.*abs(xOLOC(1,2:N+1)-xOLOC(1,1:N)));

pe=p5;

pc=pe;
%   delapx(1,:)=(piOLOC(1,:)-xOLOC(1,1:N)).*((1+t0(1,:)./xOLOC(1,1:N))); 
%   delapx(1,:)=(piOLOC(1,:)-xOLOC(1,1:N))+t0(1,:).*log(piOLOC(1,:)./xOLOC(1,1:N));  %px
%   delape(1,:)=(piOLOC(1,:)-xOLOC(1,1:N));  %pe
%   delapx(1,:)=delape(1,:).^2;
%   pc=norm(cp*miOLOC(1,:).*delat(1,:),1);
%   pc=sum(delape(1,:));
%   pc=norm((cp*miOLOC(1,:).*delapx),1)+po;
%   pc=norm(delapx(1,:),1);
%   pc=var(delat(1,:));
%   pc=sum(cp*miOLOC(1,:).*abs(delat(1,:)));
%   pc=sum((cp*miOLOC(1,:).*(piOLOC(1,N)-xOLOC(1,1:N))));
%   px=pc-sum(t0(1,:).*miOLOC(1,:)*cv*0.00001)+pm-pn+po;
%   pd=sum((cp*miOLOC(1,:).*(piOLOC(1,:)-xOLOC(1,1:N))))/COP;
%   pc=sum((1-t0(1,:)./xOLOC(1,2:N+1)).*((tdk(1,:)-xOLOC(1,2:N+1))/r12+(tdk(1,:)-xOLOC(1,2:N+1))/r13+(tdk(1,:)-xOLOC(1,2:N+1))/r14+(tdk(1,:)-xOLOC(1,2:N+1))/r15))+...
%       sum(miOLOC(1,:).*(cp*abs(piOLOC(1,:)-xOLOC(1,1:N))-cv*t0(1,:).*log(piOLOC(1,:)./xOLOC(1,1:N))))+...
%       sum(mr*(cp*(xOLOC(1,2:N+1)-xOLOC(1,1:N))-cv*t0(1,:).*log(xOLOC(1,2:N+1)./xOLOC(1,1:N))));
  
pf=sum(k*miOLOC.^3);
  
%���ʶȳͷ�
cf=(sum((pos(xOLOC(1,2:N+1)-tbk(1,:)))+(pos(tbk(2,:)-xOLOC(1,2:N+1)))))*pen;
%   cf=0;
  JOLOC = pc+pf+cf; 

   options = sdpsettings('verbose',1,'solver','fmincon');
%    options = sdpsettings();
   sol = solvesdp(F,JOLOC,options);
if sol.problem == 0
    psolution = double(piOLOC);%��ȡ�¶ȿ��Ʊ���
%     msolution = double(miOLOC);%��ȡ���ٱ���
    xsolution = double(xOLOC); %��ȡ�¶�״̬����
    fopt=double(JOLOC);        %��ȡĿ�꺯��
    
%     pe=norm((cp*miOLOC(1,:).*delape),1);
%     px=norm((cp*miOLOC(1,:).*delapx),1);
    fpt=double(pe/s); %�¶��ܺ�
    fsf=double(cf); %���ʶȳͷ�
    fpl=double(pf); %�����ܺ�
    fpx=double(px/s); %�������
%     de=double(delat);
 else
    disp('Hmm, something went wrong!');
    psolution =0;
    xsolution =0;
    yalmiperror(sol.problem)
 end
%% ==============================������=========================================
% Store the OLOC control and state vectors.
%
%   u = piOLOC;
%   x = xOLOC;

      T = psolution';
      M = miOLOC';
      x = xsolution';
      xOLOC=xsolution;
      piOLOC=psolution;
      
      T1=T-273.15;
      M1=M;
      x1=x(2:N+1,1)-273.15;
      table(T1,M1,x1)
      table(fpt,fsf,fpl,fpx)
      plot(t,tb(1,:),'cyan-.','LineWidth',2);
      hold on
      plot(t,tb(2,:),'cyan-.','LineWidth',2);
      hold on
      plot(t,x1,'b-','LineWidth',2);
      hold on
      stairs(t,T1,'r-','LineWidth',2);
      legend('���ʶ�Լ������', '���ʶ�Լ������','�����¶�', '�ͷ��¶�', 'Location', 'northeast')
      set(gca, 'XLim', [1,N], 'XTick',0:2:N,'XTickLabel',0:2:N)
      set(gca, 'YLim', [8,36], 'YTick',8:2:36,'YTickLabel',8:2:36)
      grid on
      box off
%
% �ڶ��׶�
% N =24;
% dt=1;
% t=1:dt:24;
% %% �������
% ss=[161;79;290;112;494;538;608;504;528;566;495;527;549;485;363;322;243;259;279;215;91;63;39;30]';
% tt=[22.2000   22.8000   23.7000   23.0000   22.1000 24.4571	33.3143	31.6429	28.4429	39.5143	37.1857	40.9714	42.4571	32.6714	25.0857	24.4714	25.2000   25.2000   25.1000   25.3000   25.4000   25.0000   24.8000   24.5000];
% qint=[150.8880  119.7837  118.3548   94.5784   82.8498   95.4877   87.0290  103.9159  144.2984  154.8375  300.0116  276.2081  357.5375  294.7984 354.5878  365.0463  385.5733  341.5484  325.8368  278.1317  203.6880  198.4279  182.0419  161.5226];
% td=[22.2000   22.8000   23.7000   23.0000   22.1000   23.6000   24.0000   23.2000   24.3000   25.8000   26.3000   27.5000   27.7000 26.7000   25.7000   25.7000   25.2000   25.2000   25.1000   25.3000   25.4000   25.0000   24.8000   24.5000];
% tb=[24	24	24	24	24	23	22	21	21	21	21	21	21	21	21	21	21	21  22	23	24	24	24	24;
%     17	17	17	17	17	18	19	20	20	20	20	20	20	20	20	20	20	20	19	18	17	17	17	17];
% ttk=tt+273.15;
% tdk=td+273.15;
% tbk=tb+273.15;
% % qint=0;
% %% ����¥�����
%  
%  a12=27.5;%���m2
%  a13=23;  %���m2
%  a14=27.5;%���m2
%  a15=20.2;%���m2
%  aw15=2.8;%���m2
%  
%  tw=0.9;  %����������
%  xs=0.4;   %����������
%  
%  c12=7.9e5/1000;%����kj/k
%  c13=6.6e5/1000;%����kj/k
%  c14=7.9e5/1000;%����kj/k
%  c15=2.6e7/1000;%����kj/k
%  cr=2.8e5/1000; %����kj/k
%  r12=0.0640*1000;%����k/kw
%  r13=0.0768*1000;%����k/kw
%  r14=0.0640*1000;%����k/kw
%  r15=0.0299*1000;%����k/kw
%  rw15=0.0868*1000;%����k/kw
%  
%  pen=50;   %���ʶȳͷ� kwh/k
%  pex=3.8e8;%���ʶȳͷ� kwh/k
%  h0=298.6e3/1000;%���� kj/kg
%  s0=6.86e3/1000; %���� kj/kgk
%  ft=1+273.15;    %���ʶ��ɳ�
%  
%  t0=(10+273.15)*ones(1,24);
%   
% %% �ȱò���
% md=0.52*3600;%������������kg/h
% mn=0.3*3600;    %������������kg/h
% COP=3.2;
% % m=md;
% 
% 
%  cp=1005/1000; %��ѹ������kj/kgk
%  cv=718/1000;  %�����������kj/kgk   
%  gr=287/1000;  %���峣��kj/kgk
% 
% 
% %% �������
%  k=711/(1000*3600^3);%kWh3/kg3
%  
%  ps=135;%pa
%  ad=1.29;%kg/m3
%  eff=0.15
% 
% % %% ====================================�Ż�ģ��===================================== 
% % fpc1=fpc;
% % fpf1=fpf;
% % fcf1=fcf;
% fpx1=fpx;
% %
% % Import the problem dimensions.
%  
%   nx = 5;
%   nu = 2;
%   nw = 5;
%   
%   piOLOC=sdpvar(nx,N);
% %   miOLOC=sdpvar(nx,N);
%   xOLOC =sdpvar(nx,N+1);
%   
% 
% %% ����ģ��Լ������
% 
%   % �ͷ��¶�
%   pl=16+273.15;   %�����¶�k
%   pp=16+273.15;   %�����¶�k
%   dpmax=6+273.15;
%   pmax=30+273.15;   %�¶�����k
%   pmin=17+273.15;   %�¶�����k
%   
%   %�������� 
%   ml=0.2*md;   %�����¶�k
%   mp=0.2*md;   %�����¶�k
%   dmmax=0.2*md;
%   mmax=md;    %�¶�����k
%   mmin=0;     %�¶�����k
%   
%   pMax =  [pmax; 0; 0 ;0 ;0];
%   pMin =  [pmin; 0; 0 ;0 ;0];
%   dpMax = [dpmax];
%   
%   mMax =  [mmax; 0; 0 ;0 ;0];
%   mMin =  [mmin; 0; 0 ;0 ;0];
%   dmMax = [dmmax];
%   
%   miOLOC=ones(nx,N)*md/80;
%   mr=1;
%   x0=18*ones(nx,1)+273.15;
%   qint=0;
%   %
%   % Each control vector must be in the (time-invariant) feasible control space.
%   %
%     F= (repmat(pMin,1,N) <= piOLOC <= repmat(pMax,1,N));
% %     F= F+(repmat(mMin,1,N) <= miOLOC <= repmat(mMax,1,N));
%     
%     for kp = 1 : N-1
%     F=F+  (abs(piOLOC(1,kp) - piOLOC(1,kp+1))<= dpMax);
% %     F=F+  (abs(miOLOC(1,kp) - miOLOC(1,kp+1))<= dmMax);
%     end
%   %
%   % The states, which are functions of the controls and the nominal disturbances, must obey the dynamics at all stages.
%   %
%     F=F+ (xOLOC(:,1) == x0);
%     
%     for kp = 1 : N
%     a=zeros(5,5);
%     a(1,1)=(-1/cr)*(1/r12+1/r13+1/r14+1/r15+1/rw15);
%     a(1,2)=(1/cr)*(1/r12);
%     a(1,3)=(1/cr)*(1/r13);
%     a(1,4)=(1/cr)*(1/r14);
%     a(1,5)=(1/cr)*(1/r15+1/rw15);
%     a(2,1)=(1/c12)*(1/r12);
%     a(3,1)=(1/c13)*(1/r13);
%     a(4,1)=(1/c14)*(1/r14);
%     a(5,1)=(1/c15)*(1/r15);
%     a(2,2)=(-2/c12)*(1/r12);
%     a(3,3)=(-2/c13)*(1/r13);
%     a(4,4)=(-2/c14)*(1/r14);
%     a(5,5)=(-2/c15)*(1/r15);
%     a=a+eye(5,5);
%     b=[(1/cr)*cp*miOLOC(1,kp)*(piOLOC(1,kp)-xOLOC(1,kp))];
%     w=zeros(5,24);
%     w(1,:) =(1/cr)*(tw*aw15*ss/1000+qint);
%     w(2,:) =tdk/(c12*r12);
%     w(3,:) =tdk/(c13*r13);
%     w(4,:) =tdk/(c14*r14);
%     w(5,:) =tdk/(c15*r15)+(1/c15)*xs*aw15*ss/1000;
% 
%         
%     F=F+ (xOLOC(:,kp+1) == a*xOLOC(:,kp) + b*[1;0;0;0;0] + w(:,kp));
%     end
%   %
%   % The first element of each state vector must be in the time-invariant 
%   % feasible state space. Note that the second element of the state is 
%   % unconstrained.
%   %
%    F=F+  (tbk(2,:)-1 <= xOLOC(1,2:N+1)<= tbk(1,:)+1);
%    F=F+ (piOLOC(1,:)>=xOLOC(1,1:N));
% 
%   delapx(1,:)=(piOLOC(1,:)-xOLOC(1,1:N)).*(1+t0(1,:)./xOLOC(1,1:N)); 
%   fpx=norm(cp*miOLOC(1,:).*delapx(1,:),1);
%   
%   
%   F=F+(fpx<=fpx1);
% %  F=F+(pf<=fpf1);
% %  F=F+(cf<=fcf1);
%   delat(1,:)=(piOLOC(1,:)-xOLOC(1,1:N)); 
%   pc=var(delat(1,:));
%   
%   pf=sum(k*miOLOC(1,:).^3);
%   %���ʶȳͷ�
%   cf=(sum((pos(xOLOC(1,2:N+1)-tbk(1,:)))+(pos(tbk(2,:)-xOLOC(1,2:N+1)))))*pen;
%   
% %   JOLOC=pc+cf+pf; 
%   J=pc+pf+cf;
%   
% %   JOLOC = norm(piOLOC(1,:)-xOLOC(1,2:N+1),2)+cf;
% %   JOLOC = norm((piOLOC(1,:)-xOLOC(1,1:N)),inf)+norm(((piOLOC(1,:)./t0(1,:))),inf)+cf;
% %   JOLOC = pc+norm(((piOLOC(1,:)-t0(1,:))),1)+cf;
% %   JOLOC = norm((piOLOC(1,:)+xOLOC(1,1:N))./piOLOC(1,:),1)+cf;
% %   cvx_end
%    options = sdpsettings('verbose',1,'solver','fmincon');
% %    options = sdpsettings();
%    sol = solvesdp(F,J,options);
% if sol.problem == 0
%     psolution = double(piOLOC);%��ȡ�����
%     msolution = double(miOLOC);
%     xsolution = double(xOLOC); %��ȡ�����
%     fopt=double(J); 
%     pc=norm((cp*miOLOC(1,:).*(psolution(1,:)-xsolution(1,1:N))),1);
%     fpc=double(pc);
%     fcf=double(cf);
%     fpf=double(pf);
% %     de=double(delat);
%  else
%     disp('Hmm, something went wrong!');
%     psolution =0;
%     xsolution =0;
%     yalmiperror(sol.problem)
%  end
% %% ==============================������=========================================
% % Store the OLOC control and state vectors.
% %
% %   u = piOLOC;
% %   x = xOLOC;
% 
%       T = psolution';
%       M = msolution';
%       x = xsolution';
%       xOLOC=xsolution;
%       piOLOC=psolution;
%       delat(1,:)=(psolution(1,:)-xsolution(1,1:N)).*(1+t0(1,:)./xsolution(1,1:N)); 
%       fpx=norm(cp*miOLOC(1,:).*delat(1,:),1);
% %       fpx=sum((1-t0(1,:)./xOLOC(1,2:N+1)).*((tdk(1,:)-xOLOC(1,2:N+1))/r12+(tdk(1,:)-xOLOC(1,2:N+1))/r13+(tdk(1,:)-xOLOC(1,2:N+1))/r14+(tdk(1,:)-xOLOC(1,2:N+1))/r15))+...
% %       sum(miOLOC(1,:).*(cp*(piOLOC(1,:)-xOLOC(1,1:N))-cv*t0(1,:).*log(piOLOC(1,:)./xOLOC(1,1:N))))+...
% %       sum(mr*(cp*(xOLOC(1,2:N+1)-xOLOC(1,1:N))-cv*t0(1,:).*log(xOLOC(1,2:N+1)./xOLOC(1,1:N))));
%       
%       T1=(T(:,1)-273.15);
%       M1=M(:,1);
%       x1=x(2:N+1,1)-273.15;
%       table(T1,M1,x1)
%       table(fpc,fpf,fcf,fpx)
%       figure (2)
%       plot(t,tb(1,:),'*');
%       hold on
%       plot(t,tb(2,:),'*');
%       hold on
%       plot(t,x1,'b-');
%       hold on
%       stairs(t,T1,'r-');
%       legend('���ʶ�Լ������', '���ʶ�Լ������','�����¶�', '�ͷ��¶�', 'Location', 'northeast')
%       set(gca, 'XLim', [1,N], 'XTick',0:2:N,'XTickLabel',0:2:N)
%       set(gca, 'YLim', [8,32], 'YTick',8:2:32,'YTickLabel',8:2:32)
%       grid on
%       box off
