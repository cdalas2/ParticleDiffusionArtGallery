% % % data = load("SKMC.o");
% % % N = data(:,1);
% % % t = data(:,3)/3600;
% % % 
% % % Ain = 12100.0;
% % % A = 160000.0;
% % % Aout = A - Ain;
% % % Din = 250;
% % % Dout = 3500;
% % % N_analytical = 1/(1 + Aout*Din/Ain/Dout);
% % % pe = 100*abs(N_analytical - mean(N))/N_analytical
% % % 
% % % error = 100*abs(N - N_analytical)/N_analytical;

% % % figure(3)
% % % hold on
% % % plot(t,N,'r-');
% % % plot((0:1:14500)/3600,ones(1,14501)*N_analytical,'--g');
% % % hold off
% % % 
% % % figure(4)
% % % hold on
% % % plot(t,error,'k-');
% % % hold off
% % % 
% % % figure(5)
% % % hold on
% % % for(i = 0:40)
% % %     plot(i*10*ones(1,401), 0:400, 'k-')
% % %     plot(0:400,i*10*ones(1,401), 'k-')
% % % end
% % % hold off
    

clear all
close all
clc
PositionsData = load("BilayerPDAHeadInside.o");
Sender = PositionsData(1:end-1,1)+1;
Receiver = PositionsData(1:end-1,2)+1;
rho = PositionsData(1:end-1,3);
simT = PositionsData(1:end-1,4);
K = PositionsData(end,1);
% N = ones(K,1);
N = zeros(K,1); N(Sender(1)) = 1*K;
rho_max = rho(end)
LCELLS_PER_LENGTH_SCALE = PositionsData(end,2);
TIME_MAX = PositionsData(end,4);
LENGTH_SCALE = 10; %1 micrometer
L = LENGTH_SCALE/LCELLS_PER_LENGTH_SCALE;
LatticeCoords = InitializePositionsSquare(1600,10);
Positions = zeros(N(Sender(1)),2);
% for i = 1:K
%     Positions(i,:)= LatticeCoords(i,:) + (-L/2 + (L/2 + L/2)*rand(1,2));
% end
for i = 1:N(Sender(1))
    Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,2));
end
% ParticleLocation = 1:K;
ParticleLocation = ones(N(Sender(1)),1)*Sender(1);
close all


figure(34)
hold on
% xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',19);
% ylabel('$y$ ($\mu$m)','interpreter','latex','fontsize',19);
% zlabel('$z$ ($\mu$m)','interpreter','latex','fontsize',19);
set(gca,'XColor', 'none','YColor','none')
aa = sqrt(K) * L/2;
set(gca, 'XLim', [0 2*aa], 'YLim', [0 2*aa]);
set(gca,'Color',[0.03,0.18,0.52]);
hp=plot(Positions(:,1),Positions(:,2),'.','MarkerSize',15,'MarkerEdgeColor',[0.3,0.75,0.93],...
    'MarkerFaceColor',[0.3,0.75,0.93]);
% hp=plot(Positions(:,1),Positions(:,2),'r.','MarkerSize',30);
v = VideoWriter('test');
v.Quality = 95;
open(v);
snapiter = 1;
M(snapiter) = getframe(gcf);
writeVideo(v,M(snapiter));
SNAPSHOT_RATE = 5000;
for step = 1:(length(simT)-1)
    ChosenOne = randi(N(Sender(step)));
    k = find(ParticleLocation == Sender(step),N(Sender(step)));
    ParticleLocation(k(ChosenOne)) = Receiver(step);
    N(Sender(step)) = N(Sender(step)) - 1;
    N(Receiver(step)) = N(Receiver(step)) + 1;
%     pause((simT(i+1)-simT(i))/100);
    Positions(k(ChosenOne),:) = LatticeCoords(Receiver(step),:) + (-L/2 + (L/2 + L/2)*rand(1,2));
    if mod(step,SNAPSHOT_RATE) == 0    
        set(hp,'XData',Positions(:,1),'YData',Positions(:,2));
        drawnow;
        simT(step+1)
        snapiter = snapiter + 1;
        M(snapiter) = getframe(gcf);
        writeVideo(v,M(snapiter));
    end
end

close(v);
% movie(M);
hold off

