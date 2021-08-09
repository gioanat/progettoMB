% clc
% clear
% load datiFittingCovid19.mat
% 
% tStart = datetime(2020,11,04);
% tEnd = datetime(2020,11,14);
% indT = find(time>=tStart & time <=tEnd);
% dt = 1/24; % time step
% temp=time(indT);
% timeA=datetime(temp(1), 'Locale', 'en_US'):dt:datetime(temp(end), 'Locale', 'en_US');
% N = numel(timeA);
% t = [0:N-1].*dt;
% 
% A=11*0.85*[0 118e03/popolazione(2) 0 104e03/popolazione(4) 0 0;
%     118e03/popolazione(1) 0 286e03/popolazione(3) 283e03/popolazione(4) 0 0;
%     0 286e03/popolazione(2) 0 129e03/popolazione(4) 0 0;
%     104e03/popolazione(1) 283e03/popolazione(2) 129e03/popolazione(3) 0 286e03/popolazione(5)  109e03/popolazione(6);
%     0 0 0 286e03/popolazione(4) 0 0;
%     0 0 0 109e03/popolazione(4) 0 0];
% 
% flag=1; % diffusion yes
% 
% indexPeriod=6;
% 
% S0 = squeeze(SEIQRDP{indexPeriod-1}.YM(1,:,end));
% P0 = squeeze(SEIQRDP{indexPeriod-1}.YM(2,:,end));
% E0 = squeeze(SEIQRDP{indexPeriod-1}.YM(3,:,end));
% I0 = squeeze(SEIQRDP{indexPeriod-1}.YM(4,:,end));
% Q0 = squeeze(SEIQRDP{indexPeriod-1}.YM(5,:,end));
% R0 = squeeze(SEIQRDP{indexPeriod-1}.YM(6,:,end));
% D0 = squeeze(SEIQRDP{indexPeriod-1}.YM(7,:,end));
% 
% epsilon = [0.01:0.01:0.1 0.1:0.1:1];
% 
% num_edge = 14;
% top = dec2bin(2^num_edge-1:-1:0);
% 
% for i = 1:num_edge
%     topologie(:,i) =str2num(top(:,i)); 
% end
% 
% for i=1:2^num_edge
%     if mod(i,100) == 0
%         disp(int2str(i/(2^num_edge)*100))
%     end
%     M = zeros(6,6);
%     M(find(A>0)) = topologie(i,:)'.*A(find(A>0));
%     Ldiag=diag(sum(M,1));
%     L=M-Ldiag;
%     for j=1:length(epsilon)
%         Y = simulatedSEIQRDP_onNetwork(popolazione,S0,P0,E0,I0,Q0,R0,D0,Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),L,[SEIQRDP_onNetwork{indexPeriod-1}.parametriM(1:end-2) epsilon(j) epsilon(j)],N,t,flag);
%         mat(i,j) = mean(sum(Y(5,:,:)));
%     end
% end

clc
clear
close all
load topologie.mat
figure
hold on
plot(1:2^(num_edge),mat(:,end),'.',1:2^(num_edge),mat(end,end)*ones(size(mat(:,end))),'k-','LineWidth',2)
[mm, idx] = max(mat(:,end));
plot(idx,mm,'ro')
[mm, idx] = min(mat(:,end));
plot(idx,mm,'go')
hold off
xlabel("Topologia")
ylabel("Numero medio giornaliero stimato della categoria Q")
axis([1 2^(14) min(mat(:,end))-1000 max(mat(:,end))+1000])
 
figure
% filename = 'testnew51.gif';
s = [1 1 2 2 2 3 3 4 4 4 4 4 5 6];
t = [2 4 1 3 4 2 4 1 2 3 5 6 4 4];
for i= [1 52 1013 497 13197 3447]
      subplot(1,2,1)
      plot(epsilon,mat(end,:),epsilon,mat(i,:))
      legend("Topologia vuota",strcat("Topologia ",int2str(i)))
      title(strcat("Topologia ",int2str(i)))
      xlabel("\epsilon")
      ylabel("Numero medio giornaliero stimato della categoria Q")
      subplot(1,2,2)
      temp = topologie(i,:)'.*A(find(A>0));
      weight = temp(logical(topologie(i,:)));
      G = digraph(s(logical(topologie(i,:))),t(logical(topologie(i,:))),weight);
      plot(G)
      title(strcat("Topologia ",int2str(i)))
      pause
      drawnow
      frame = getframe(1);
      im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if n == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%       end
end
epsilon = [0.01:0.01:0.1 0.1:0.1:1];
M=11*[0 0.55*118e03/popolazione(2) 0 0.55*104e03/popolazione(4) 0 0;
    0.55*118e03/popolazione(1) 0 0.55*286e03/popolazione(3) 0.55*283e03/popolazione(4) 0 0;
    0 0.55*286e03/popolazione(2) 0 0.85*129e03/popolazione(4) 0 0;
    0.55*104e03/popolazione(1) 0.55*283e03/popolazione(2) 0.85*129e03/popolazione(3) 0 0.85*286e03/popolazione(5)  0.85*109e03/popolazione(6);
    0 0 0 0.85*286e03/popolazione(4) 0 0;
    0 0 0 0.85*109e03/popolazione(4) 0 0];
Ldiag=diag(sum(M,1));
L=M-Ldiag;
for j=1:length(epsilon)
    Y = simulatedSEIQRDP_onNetwork(popolazione,S0,P0,E0,I0,Q0,R0,D0,Qmedia(indT,:),Rmedia(indT,:),Dmedia(indT,:),L,[SEIQRDP_onNetwork{indexPeriod-1}.parametriM(1:end-2) epsilon(j) epsilon(j)],N,t,flag);
    vigente(j) = mean(sum(Y(5,:,:)));
end
figure
subplot(1,2,1)
plot(epsilon,mat(end,:),epsilon,vigente)
legend("Topologia vuota","Topologia (politica vigente)")
title("Topologia (politica vigente)")
xlabel("\epsilon")
ylabel("Numero medio giornaliero stimato della categoria Q")
subplot(1,2,2)
temp = topologie(1,:)'.*M(find(M>0));
weight = temp(logical(topologie(1,:)));
eLabels = {'-30%' '-30%' '-30%' '-30%' '-30%' '-30%' '' '-30%' '-30%' '' '' '' '' ''};
G = digraph(s(logical(topologie(1,:))),t(logical(topologie(1,:))),weight);
nodeColor = [1 0 0;
    1 0 0;
    1 1 0;
    1 1 0;
    1 1 0;
    1 1 0];
plot(G,'NodeColor',nodeColor,'EdgeLabel',eLabels)

title(strcat("Topologia ",int2str(1)))
