%%plot IVIV
% load IVIVcompare
clear IVIV_Coms IVIV_TLM IVIV_Comsol ID t IVIV_Ansys
Num=501;
ID=1:Num;
Delt=0.1;
t=ID*Delt;

IVIV_TLM=dlmread('IVIV.txt','',0,0);
%IVIV_TLM1=dlmread('IVIV_Comsol1.txt','',8,0);
IVIV_Comsol=dlmread('IVIV_Comsol.txt','',8,0);

% IVIV_TLM=dlmread('D:\Backup\Peng\Peng4thPaper\SinglePhase Transformer\SinglePhaseTransformer\IVIV.txt','',0,0);
% IVIV_Comsol=dlmread('D:\Backup\Peng\Peng4thPaper\SinglePhase Transformer\Comsol\IVIV_Comsol.txt','',0,0);

IVIV_Coms(:,1)=IVIV_Comsol(1:Num,2);
IVIV_Coms(:,3)=IVIV_Comsol(Num+1:Num*2,2);
IVIV_Coms(:,2)=IVIV_Comsol(Num*2+1:Num*3,2);
IVIV_Coms(:,4)=IVIV_Comsol(Num*3+1:Num*4,2);

IVIV_TLM(:,1)=IVIV_TLM1(1:Num,2);
IVIV_TLM(:,3)=IVIV_TLM1(Num+1:Num*2,2);
IVIV_TLM(:,2)=IVIV_TLM1(Num*2+1:Num*3,2);
IVIV_TLM(:,4)=IVIV_TLM1(Num*3+1:Num*4,2);


figure(1)
set(gcf,'pos',[0 0 800 500])
subplot_tight(4,1,1);
hold on;
plot(t,1.414*IVIV_TLM(ID,1)','r--','LineWidth',2);
plot(t,1.414*IVIV_Coms(ID,1)','b-','LineWidth',1);
% xlabel('Time (ms)','FontSize',8.5,'FontName','Arial');
ylabel('I_p (A)','FontSize',8.5,'FontName','Arial');
set(gca,'FontSize',8.5,'FontName','Arial');
legend('NDD','Comsol');
legend boxoff
axis([0 100 -2000 2000]);
box on;

subplot_tight(4,1,2);
hold on;
plot(t,1.414*IVIV_TLM(ID,3)','r--','LineWidth',2);
plot(t,1.414*IVIV_Coms(ID,3)','b-','LineWidth',1);
% legend('NDD','Comsol');
% xlabel('Time (ms)','FontSize',8.5,'FontName','Arial');
ylabel('I_s (A)','FontSize',8.5,'FontName','Arial');
set(gca,'FontSize',8.5,'FontName','Arial');
axis([0 100 -6000 6000]);
% axis auto;
box on;
subplot_tight(4,1,3);
hold on;
plot(t,1.414*IVIV_TLM(ID,2)/1000','r--','LineWidth',2);
plot(t,1.414*IVIV_Coms(ID,2)/1000','b-','LineWidth',1);
% xlabel('Time (ms)','FontSize',8.5,'FontName','Arial');
ylabel('U_p (kV)','FontSize',8.5,'FontName','Arial');
set(gca,'FontSize',8.5,'FontName','Arial');
% legend('NDD','Comsol');
axis([0 100 -1.414*2.5e2 1.414*2.5e2]);
% axis auto;
box on;

subplot_tight(4,1,4);
hold on;
plot(t,1.414*IVIV_TLM(ID,4)/1000','r--','LineWidth',2);
plot(t,1.414*IVIV_Coms(ID,4)/1000','b-','LineWidth',1);
% legend('NDD','Comsol');
xlabel('Time (ms)','FontSize',8.5,'FontName','Arial');
ylabel('U_s (kV)','FontSize',8.5,'FontName','Arial');
set(gca,'FontSize',8.5,'FontName','Arial');
axis([0 100 -4000 4000]);
% axis auto;
box on;





