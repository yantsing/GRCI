close all

linewidth = 2
frontsize =14

load Phy

figure(1)
box on
hold on
plot(Cs, squeeze(ovals_2(1,:)), '-*b','LineWidth',linewidth)
hold on
plot(Cs, squeeze(ovals_3(1,:)), '--or','LineWidth',linewidth)
plot(Cs, squeeze(ovals_4(1,:)), '-.^g','LineWidth',linewidth)

legend('N = 2','N = 3','N = 4')
xlabel('Worst-case secrecy rate R_s (b/s/Hz)','fontsize',frontsize)
ylabel('Total transmit power','fontsize',frontsize)
grid on
pause
print2eps Changing_Cs
eps2pdf Changing_Cs.eps Changing_Cs.pdf


load Phy

figure(2)
plot(Ss, squeeze(ovals_2(:,6)), '-*b','LineWidth',linewidth)
hold on
plot(Ss, squeeze(ovals_3(:,6)), '--or','LineWidth',linewidth)
plot(Ss, squeeze(ovals_4(:,6)), '-.^g','LineWidth',linewidth)
set(gca,'XLim',[1 20])
set(gca,'XTick',[1,5,10,15,20])

legend('N = 2','N = 3','N = 4')
xlabel('Simulation parameter s','fontsize',frontsize)
ylabel('Total transmit power','fontsize',frontsize)
grid on
pause
print2eps Changing_s
eps2pdf Changing_s.eps Changing_s.pdf



Cs_Success= [success2(1,:);success3(1,:);success4(1,:)];
Cs_Success = Cs_Success';
Cs_Success = Cs_Success/ max_trial;


load lesshelper
figure(3)
% plot(Cs,Cs_Success(:,1),'-*b','LineWidth',linewidth);
% hold on
% plot(Cs,Cs_Success(:,2),'--or','LineWidth',linewidth);
% plot(Cs,Cs_Success(:,3),'-.^g','LineWidth',linewidth);

plot(Cs,success2(1,:),'-*b','LineWidth',linewidth);
hold on
plot(Cs,success3(1,:),'--or','LineWidth',linewidth);
plot(Cs,success4(1,:),'-.^g','LineWidth',linewidth);

legend('N = 2','N = 3','N = 4')
xlabel('Worst-case secrecy rate R_s (b/s/Hz)','fontsize',frontsize)
ylabel({'Probability of satisfying'; 'the secrecy rate requirement'},'fontsize',frontsize)
axis([0.5 4 0.7 1]) 
yt = get(gca, 'ytick'); ytl = strcat(strtrim(cellstr(num2str(yt'*100))), '%'); set(gca, 'yticklabel', ytl); 
grid on
pause
print2eps Cs_success
eps2pdf Cs_success.eps Cs_success.pdf


s_Success= [success2(:,6),success3(:,6),success4(:,6)];
s_Success = s_Success/ max_trial;

figure(4)
load lesshelper
plot(Ss,success2(:,6),'-*b','LineWidth',linewidth);
hold on
plot(Ss,success3(:,6),'--or','LineWidth',linewidth);
plot(Ss,success4(:,6),'-.^g','LineWidth',linewidth);
legend('N = 2','N = 3','N = 4')
xlabel('Simulation paramenter s','fontsize',frontsize)
ylabel({'Probability of satisfying'; 'the secrecy rate requirement'},'fontsize',frontsize)
axis([1 20 0.7 1])
yt = get(gca, 'ytick'); ytl = strcat(strtrim(cellstr(num2str(yt'*100))), '%'); set(gca, 'yticklabel', ytl); 
grid on
pause
print2eps s_success
eps2pdf s_success.eps s_success.pdf


% figure(5)
% [m,n] = size(Cs_Success);
% Cs_outage = zeros(m,n);
% for i = 1: m
%     for j = 1: n
%         Cs_outage(i,j) = 1- Cs_Success(i,j);
%     end
% end
% plot(Cs,Cs_outage(:,1),'-*b','LineWidth',linewidth);
% hold on
% plot(Cs,Cs_outage(:,2),'--or','LineWidth',linewidth);
% plot(Cs,Cs_outage(:,3),'-.^g','LineWidth',linewidth);
% 
% legend('N = 2','N = 3','N = 4')
% xlabel('Required secure rate R_s (b/s/Hz)','fontsize',frontsize)
% ylabel('Secure communications outage rate','fontsize',frontsize)
% axis([0.5 4 0.3 0]) 
% yt = get(gca, 'ytick'); ytl = strcat(strtrim(cellstr(num2str(yt'*100))), '%'); set(gca, 'yticklabel', ytl); 
% grid on
% pause
% print2eps Cs_outage
% eps2pdf Cs_outage.eps Cs_outage.pdf



% figure(6)
% 
% [m,n] = size(s_Success);
% for i = 1 : m
%     for j = 1 : n
%         s_outage = 1- s_Success(i,j);
%     end
% end
% plot(Ss,s_outage(:,1),'-*b','LineWidth',linewidth);
% hold on
% plot(Ss,s_outage(:,2),'--or','LineWidth',linewidth);
% plot(Ss,s_outage(:,3),'-.^g','LineWidth',linewidth);
% legend('N = 2','N = 3','N = 4')
% xlabel('Simulation paramenter s','fontsize',frontsize)
% ylabel('Secure communications outage rate','fontsize',frontsize)
% axis([1 20 0.3 0])
% yt = get(gca, 'ytick'); ytl = strcat(strtrim(cellstr(num2str(yt'*100))), '%'); set(gca, 'yticklabel', ytl); 
% grid on
% pause
% print2eps s_outage
% eps2pdf s_outage.eps s_outage.pdf


figure(16)
load lesshelper
box on
hold on
plot(Cs, success4(1,:), '-bv','LineWidth',linewidth)
plot(Cs, success4(2,:), '-b^','LineWidth',linewidth)
plot(Cs, success0(1,:), '--kv','LineWidth',linewidth)
plot(Cs, success0(2,:), '--k^','LineWidth',linewidth)
plot(Cs, success1(1,:), '-.rv','LineWidth',linewidth)
plot(Cs, success1(2,:), '-.r^','LineWidth',linewidth)


legend('s = 1','s = 5')
xlabel('Worst-case secrecy rate R_s (b/s/Hz)','fontsize',frontsize)
ylabel({'Probability of satisfying'; 'the secrecy rate requirement'},'fontsize',frontsize)
% axis([1 20 0.7 1])
yt = get(gca, 'ytick'); ytl = strcat(strtrim(cellstr(num2str(yt'*100))), '%'); set(gca, 'yticklabel', ytl); 
grid on
pause
print2eps comparison
eps2pdf comparison.eps comparison.pdf



