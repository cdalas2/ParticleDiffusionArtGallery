% data = load("SKMCTM3D_2.o");
% pt = data(1:(end-1),1);
% rho = data(1:(end-1),2);
% Sender = data(1:(end-1),3);
% Receiver = data(1:(end-1),4);
% t = data(1:(end-1),5)/3600;
% 
% SNAPSHOT_RATE = data(end,1);
% SDSPEED = data(end,2);
% TOTAL_LATTICE_CELLS = data(end,3);
% SCALE = data(end,4);
% TIME_MAX = data(end,5);
% 
Vin = (180^2)*400;
V = 400^3;
Vout = V - Vin;
Din = 180;
Dout = 540;
% rhomean = mean(rho((end-100):end))
% rhoav = rhomean*ones(14501,1);
rho_upperbound = ones(14501,1)/(1 + Vout*Din/Vin/Dout)
% 
% 
figure(3)
hold on
% plot(t,rho,'r-');
% plot((0:1:14500)/3600,rhoav,'--g');
plot((0:1:14500)/3600,rho_upperbound,'k.');
hold off

% % % clear all
% % % close all
% % % clc
% % % PositionsData = load("SKMCTM3D_MovieData.o");
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
% % % LatticeCoords = InitializePositionsCube(K,L);
% % % 
% % % Positions = LatticeCoords;
% % % % for i = 1:N(Sender(1))
% % % %     Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
% % % % end
% % % ParticleLocation = 1:K;
% % % SDBase = 1;
% % % DL = K^(1/3);
% % % MDL = 18;
% % % % ParticleLocation = ones(N(Sender(1)),1)*Sender(1);
% % % 
% % % figure(34)
% % % hold on
% % % xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',19);
% % % ylabel('$y$ ($\mu$m)','interpreter','latex','fontsize',19);
% % % zlabel('$z$ ($\mu$m)','interpreter','latex','fontsize',19);
% % % 
% % % view(43,24);
% % % aa = DL * L/2;
% % % set(gca, 'XLim', [-aa aa], 'YLim', [-aa aa], 'ZLim', [-aa aa]);
% % % set(gca,'Color','none');
% % % hp=plot3(Positions(:,1),Positions(:,2),Positions(:,3),'b.','MarkerSize',30);
% % % v = VideoWriter('BacteriaBusRide');
% % % v.Quality = 95;
% % % open(v);
% % % M(1) = getframe(gcf);
% % % writeVideo(v,M(1));
% % % for step = 1:(length(simT)-1)
% % %     if Sender(step) ~= Receiver(step)
% % %         ChosenOne = randi(N(Sender(step)));
% % %         k = find(ParticleLocation == Sender(step),N(Sender(step)));
% % %         ParticleLocation(k(ChosenOne)) = Receiver(step);
% % %         N(Sender(step)) = N(Sender(step)) - 1;
% % %         N(Receiver(step)) = N(Receiver(step)) + 1;
% % % %     pause((simT(i+1)-simT(i))/100);
% % %         Positions(k(ChosenOne),:) = LatticeCoords(ParticleLocation(k(ChosenOne)),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
% % %         set(hp,'XData',Positions(:,1),'YData',Positions(:,2),'ZData',Positions(:,3));
% % %     else
% % %         for k = 1:MDL
% % %             for i = 1:MDL
% % %                 for j = (MDL+SDBase):(1+SDBase)
% % %                     newj = j;
% % %                     newj1 = j-1;
% % %                     if j>DL
% % %                         newj = j-DL;
% % %                     end
% % %                     if (j-1) > DL
% % %                         newj1 = j-1-DL;
% % %                     end
% % %                     if j ~=(MDL+SDBase)
% % %                         N(k*DL^2 + i*DL + newj) = N(k*DL^2 + i*DL + newj1);
% % %                     else
% % %                         N(k*DL^2 + i*DL + newj) = N(k*DL^2 + i*DL + newj) + N(k*DL^2 + i*DL + newj1);
% % %                     end
% % %                     N(k*DL^2 + i*DL + SDBase) = 0;
% % %                     ParticleLocation(k*DL^2 + i*DL + newj1) = ParticleLocation(k*DL^2 + i*DL + newj);
% % %                     Positions(k*DL^2 + i*DL + newj1,:) = LatticeCoords(ParticleLocation(k*DL^2 + i*DL + newj1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
% % %                     set(hp,'XData',Positions(:,1),'YData',Positions(:,2),'ZData',Positions(:,3));
% % %                 end
% % %             end
% % %         end
% % %         SDBase = SDBase + 1;
% % %     end
% % %     drawnow;
% % %     simT(step+1)
% % %     M(step+1) = getframe(gcf);
% % %     writeVideo(v,M(step+1));
% % %     view(43 + 0.02*step,24+ 0.02*step);
% % % end
% % % close(v);
% % % movie(M);
% % % hold off
% % % 
