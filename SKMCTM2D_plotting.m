% % % data = load("SKMCTM2D.o");
% % % pt = data(1:(end-1),1);
% % % rho = data(1:(end-1),2);
% % % Sender = data(1:(end-1),3);
% % % Receiver = data(1:(end-1),4);
% % % t = data(1:(end-1),5)/3600;
% % % 
% % % SNAPSHOT_RATE = data(end,1);
% % % SDSPEED = data(end,2);
% % % TOTAL_LATTICE_CELLS = data(end,3);
% % % SCALE = data(end,4);
% % % TIME_MAX = data(end,5);

clear all
clc
SquarePositionsData = load("SKMC2.o");
SquareSender = SquarePositionsData(1:end-1,1)+1;
SquareReceiver = SquarePositionsData(1:end-1,2)+1;
Squarerho = SquarePositionsData(1:end-1,3);
SquaresimT = SquarePositionsData(1:end-1,4);
SquareK = SquarePositionsData(end,1);
SquareN = ones(K,1);

Squarerho_max = Squarerho(end)
SquareLCELLS_PER_LENGTH_SCALE = SquarePositionsData(end,2);
SquareTIME_MAX = SquarePositionsData(end,4);
SquareLENGTH_SCALE = 10; %1 micrometer
SquareL = SquareLENGTH_SCALE/SquareLCELLS_PER_LENGTH_SCALE;
SquareLatticeCoords = InitializePositionsSquare(1600,10);
for i = 1:SquareK
    SquarePositions(i,:)= SquareLatticeCoords(i,:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
end
SquareParticleLocation = 1:SquareK;

Ain = 180*400;
A = 400^2;
Aout = A - Ain;
Din = 180;
Dout = 540;
rhomean = mean(rho((end-100):end))
rhoav = rhomean*ones(14501,1);
rho_upperbound = ones(14501,1)/(1 + Aout*Din/Ain/Dout)


figure(3)
hold on
% plot(t,rho,'r-');
% plot((0:1:14500)/3600,rhoav,'--g');
plot((0:1:14500)/3600,rho_upperbound,'k.');
hold off

% % % clear all
% % % close all
% % % clc
% % % PositionsData = load("SKMCTM2D_MovieData.o");
% % % Sender = PositionsData(1:end-1,1)+1;
% % % Receiver = PositionsData(1:end-1,2)+1;
% % % simT = PositionsData(1:end-1,3);
% % % K = PositionsData(end,1);
% % % 
% % % N = ones(K,1);
% % % % N = zeros(K,1); N(Sender(1)) = 200;
% % % 
% % % 
% % % SDSPEED = PositionsData(end,2);
% % % TIME_MAX = PositionsData(end,3);
% % % LCELLS_PER_LENGTH_SCALE = 1;
% % % LENGTH_SCALE = 10; %10 micrometers
% % % L = LENGTH_SCALE/LCELLS_PER_LENGTH_SCALE;
% % % LatticeCoords = InitializePositionsSquare(K,L);
% % % 
% % % Positions = LatticeCoords;
% % % % for i = 1:N(Sender(1))
% % % %     Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
% % % % end
% % % ParticleLocation = (1:K)';
% % % SDBase = 1;
% % % DL = K^(1/2);
% % % MDL = 18;
% % % % ParticleLocation = ones(N(Sender(1)),1)*Sender(1);
% % % 
% % % figure(34)
% % % hold on
% % % xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',19);
% % % ylabel('$y$ ($\mu$m)','interpreter','latex','fontsize',19);
% % % 
% % % aa = DL * L/2;
% % % set(gca, 'XLim', [0 2*aa], 'YLim', [0 2*aa]);
% % % set(gca,'Color','none');
% % % hp=plot(Positions(:,1),Positions(:,2),'b.','MarkerSize',30);
% % % v = VideoWriter('BacteriaBusRide');
% % % v.Quality = 95;
% % % open(v);
% % % M(1) = getframe(gcf);
% % % writeVideo(v,M(1));
% % % for step = 1:(length(simT))
% % %     if Sender(step) ~= Receiver(step)
% % %         ChosenOne = randi(N(Sender(step)));
% % %         k = find(ParticleLocation == Sender(step),N(Sender(step)));
% % %         ParticleLocation(k(ChosenOne)) = Receiver(step);
% % %         N(Sender(step)) = N(Sender(step)) - 1;
% % %         N(Receiver(step)) = N(Receiver(step)) + 1;
% % % %     pause((simT(i+1)-simT(i))/100);
% % %         Positions(k(ChosenOne),:) = LatticeCoords(ParticleLocation(k(ChosenOne)),:) + (-L/2 + (L/2 + L/2)*rand(1,2));
% % %         set(hp,'XData',Positions(:,1),'YData',Positions(:,2));
% % %     else       
% % %         yell = 1
% % %         for i = 0:(MDL-1)
% % %             for j = (MDL+SDBase):(-1):(1+SDBase)
% % %                 newj = j + i*DL;
% % %                 newj1 = j-1 + i*DL;
% % %                 if j>DL
% % %                     newj = j-DL + i*DL;
% % %                 end
% % %                 if (j-1) > DL
% % %                     newj1 = j-1-DL + i*DL;
% % %                 end
% % %                 if N(newj1) ~= 0
% % %                     k = find(ParticleLocation == newj1,N(newj1));
% % %                     for p = 1:length(k)
% % %                         ParticleLocation(k(p)) = newj;
% % %                         N(newj1) = N(newj1) - 1;
% % %                         N(newj) = N(newj) + 1;
% % %                         Positions(k(p),:) = LatticeCoords(newj,:) + (-L/2 + (L/2 + L/2)*rand(1,2));
% % %                     end
% % %                 end
% % % %                 if j ~= (MDL+SDBase)
% % % %                     N(i*DL + newj) = N(i*DL + newj1);
% % % %                 else
% % % %                     N(i*DL + newj) = N(i*DL + newj) + N(i*DL + newj1);
% % % %                 end
% % %             end
% % % %             N(i*DL + SDBase) = 0;
% % %         end
% % %         if SDBase <= DL
% % %             SDBase = SDBase + 1;
% % %         else
% % %             SDBase = SDBase + 1 - DL;
% % %         end
% % %         set(hp,'XData',Positions(:,1),'YData',Positions(:,2));
% % %     end
% % %     if mod(step,333) == 0
% % %         drawnow;
% % %         simT(step)
% % %         M(step+1) = getframe(gcf);
% % %         writeVideo(v,M(step+1));
% % %     end
% % % end
% % % close(v);
% % % movie(M);
% % % hold off
% % % 
