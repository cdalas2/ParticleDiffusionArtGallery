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
close all
clc
DIM = 2;
SquarePositionsData = load("BilayerPDAHeadOutline.o");
SquareSender = SquarePositionsData(1:end-1,1)+1;
SquareReceiver = SquarePositionsData(1:end-1,2)+1;
Squarerho = SquarePositionsData(1:end-1,3);
SquaresimT = SquarePositionsData(1:end-1,4);
SquareK = SquarePositionsData(end,1);

% SquareN = ones(SquareK,1);
SquareN = zeros(SquareK,1); SquareN(SquareSender(1)) = SquareK;
Squarerho_max = Squarerho(end)
SquareLCELLS_PER_LENGTH_SCALE = SquarePositionsData(end,2);
SquareTIME_MAX = SquarePositionsData(end,4);
SquareLENGTH_SCALE = 10; %1 micrometer
SquareL = SquareLENGTH_SCALE/SquareLCELLS_PER_LENGTH_SCALE;
SquareLatticeCoords = InitializePositionsSquare(1600,10);
SquarePositions = zeros(SquareN(SquareSender(1)),DIM);
for i = 1:SquareN(SquareSender(1))
%     SquarePositions(i,:)= SquareLatticeCoords(i,:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
      SquarePositions(i,:)= SquareLatticeCoords(SquareSender(1),:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
end
% SquareParticleLocation = 1:SquareK;
SquareParticleLocation = ones(SquareN(SquareSender(1)),1)*SquareSender(1);

% % % Ain = 180*400;
% % % A = 400^2;
% % % Aout = A - Ain;
% % % Din = 180;
% % % Dout = 540;
% % % rhomean = mean(rho((end-100):end))
% % % rhoav = rhomean*ones(14501,1);
% % % rho_upperbound = ones(14501,1)/(1 + Aout*Din/Ain/Dout)


% % % figure(3)
% % % hold on
% % % % plot(t,rho,'r-');
% % % % plot((0:1:14500)/3600,rhoav,'--g');
% % % plot((0:1:14500)/3600,rho_upperbound,'k.');
% % % hold off

PositionsData = load("BilayerPDAHeadInside.o");
Sender = PositionsData(1:end-1,1)+1;
Receiver = PositionsData(1:end-1,2)+1;
rho = PositionsData(1:end-1,3);
simT = PositionsData(1:end-1,4);
K = PositionsData(end,1);
% N = ones(K,1);
N = zeros(K,1); N(Sender(1)) = K;
rho_max = rho(end)
LCELLS_PER_LENGTH_SCALE = PositionsData(end,2);
TIME_MAX = PositionsData(end,4);
LENGTH_SCALE = 10; %1 micrometer
L = LENGTH_SCALE/LCELLS_PER_LENGTH_SCALE;
LatticeCoords = InitializePositionsSquare(1600,10);
Positions = zeros(N(Sender(1)),DIM);
for i = 1:N(Sender(1))
%     Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
    Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,2));

end
% ParticleLocation = (1:K)';
SDBase = 1;
DL = K^(1/2);
MDL = 18;
ParticleLocation = ones(N(Sender(1)),1)*Sender(1);

figure(34)
hold on
xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',19);
ylabel('$y$ ($\mu$m)','interpreter','latex','fontsize',19);
set(gca,'XColor', 'none','YColor','none')
aa = DL * L/2;
set(gca, 'XLim', [0 2*aa], 'YLim', [0 2*aa]);
set(gca,'Color','none');
% set(gca,'Color',[153 27 30]/255);
% hp=plot(Positions(:,1),Positions(:,2),'.','MarkerSize',30,'MarkerEdgeColor',[255,204,0]/255,...
%     'MarkerFaceColor',[255,204,0]/255);
hp=plot(Positions(:,1),Positions(:,2),'b.','MarkerSize',20);
Squarehp=plot(SquarePositions(:,1),SquarePositions(:,2),'r.','MarkerSize',20);
% comphp=plot([SquarePositions(:,1);Positions(:,1)],[SquarePositions(:,2);Positions(:,2)],'r.','MarkerSize',30);
v = VideoWriter('BilayerPDA');
v.Quality = 95;
open(v);
M(1) = getframe(gcf);
writeVideo(v,M(1));

if (length(SquaresimT) >= length(simT))
    Tsteps = length(simT);
elseif (length(SquaresimT) < length(simT))
    Tsteps = length(SquaresimT);
end

for step = 1:Tsteps
    
    SquareChosenOne = randi(SquareN(SquareSender(step)));
    Squarek = find(SquareParticleLocation == SquareSender(step),SquareN(SquareSender(step)));
    SquareParticleLocation(Squarek(SquareChosenOne)) = SquareReceiver(step);
    SquareN(SquareSender(step)) = SquareN(SquareSender(step)) - 1;
    SquareN(SquareReceiver(step)) = SquareN(SquareReceiver(step)) + 1;
%     pause((simT(i+1)-simT(i))/100);
    SquarePositions(Squarek(SquareChosenOne),:) = SquareLatticeCoords(SquareReceiver(step),:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
    
        ChosenOne = randi(N(Sender(step)));
        k = find(ParticleLocation == Sender(step),N(Sender(step)));
        ParticleLocation(k(ChosenOne)) = Receiver(step);
        N(Sender(step)) = N(Sender(step)) - 1;
        N(Receiver(step)) = N(Receiver(step)) + 1;
%     pause((simT(i+1)-simT(i))/100);
        Positions(k(ChosenOne),:) = LatticeCoords(ParticleLocation(k(ChosenOne)),:) + (-L/2 + (L/2 + L/2)*rand(1,2));
   set(hp,'XData',Positions(:,1),'YData',Positions(:,2));
    if mod(step,5000) == 0  
        drawnow;
    end
%     set(comphp,'XData',[SquarePositions(:,1);Positions(:,1)],'YData',[SquarePositions(:,2);Positions(:,2)]);
    set(Squarehp,'XData',SquarePositions(:,1),'YData',SquarePositions(:,2));    
    if mod(step,100) == 0
        drawnow;
        simT(step)
        M(step+1) = getframe(gcf);
        writeVideo(v,M(step+1));
    end
end
% drawnow
close(v);
movie(M);
hold off

% close(v);
% % movie(M);
% hold off

