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
SquarePositionsData = load("Comp1.o");
SquareSender = SquarePositionsData(1:end-1,1)+1;
SquareReceiver = SquarePositionsData(1:end-1,2)+1;
Squarerho = SquarePositionsData(1:end-1,3);
SquaresimT = SquarePositionsData(1:end-1,4);
SquareK = SquarePositionsData(end,1);
SquareN = ones(SquareK,1);
% SquareN = zeros(SquareK,1); SquareN(SquareSender(1)) = SquareK;
Squarerho_max = Squarerho(end)
SquareLCELLS_PER_LENGTH_SCALE = SquarePositionsData(end,2);
SquareTIME_MAX = SquarePositionsData(end,4);
SquareLENGTH_SCALE = 10; %1 micrometer
SquareL = SquareLENGTH_SCALE/SquareLCELLS_PER_LENGTH_SCALE;
SquareLatticeCoords = InitializePositionsSquare(1600,10);
SquarePositions = zeros(size(SquareLatticeCoords));
for i = 1:SquareK
    SquarePositions(i,:)= SquareLatticeCoords(i,:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
%       SquarePositions(i,:)= SquareLatticeCoords(SquareSender(1),:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
end
SquareParticleLocation = 1:SquareK;
% SquareParticleLocation = ones(SquareN(SquareSender(1)),1)*SquareSender(1);

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

PositionsData = load("SKMCTM2D_MovieData.o");
Sender = PositionsData(1:end-1,1)+1;
Receiver = PositionsData(1:end-1,2)+1;
simT = PositionsData(1:end-1,3);
K = PositionsData(end,1);

N = ones(K,1);
% N = zeros(K,1); N(Sender(1)) = 200;


SDSPEED = PositionsData(end,2);
TIME_MAX = PositionsData(end,3);
LCELLS_PER_LENGTH_SCALE = 1;
LENGTH_SCALE = 10; %10 micrometers
L = LENGTH_SCALE/LCELLS_PER_LENGTH_SCALE;
LatticeCoords = InitializePositionsSquare(K,L);

Positions = LatticeCoords;
% for i = 1:N(Sender(1))
%     Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
% end
ParticleLocation = (1:K)';
SDBase = 1;
DL = K^(1/2);
MDL = 18;
% ParticleLocation = ones(N(Sender(1)),1)*Sender(1);

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
hp=plot(Positions(:,1),Positions(:,2),'b.','MarkerSize',30);
Squarehp=plot(SquarePositions(:,1),SquarePositions(:,2),'r.','MarkerSize',30);
% comphp=plot([SquarePositions(:,1);Positions(:,1)],[SquarePositions(:,2);Positions(:,2)],'r.','MarkerSize',30);
v = VideoWriter('CompositionPiece1');
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
    
    if Sender(step) ~= Receiver(step)
        ChosenOne = randi(N(Sender(step)));
        k = find(ParticleLocation == Sender(step),N(Sender(step)));
        ParticleLocation(k(ChosenOne)) = Receiver(step);
        N(Sender(step)) = N(Sender(step)) - 1;
        N(Receiver(step)) = N(Receiver(step)) + 1;
%     pause((simT(i+1)-simT(i))/100);
        Positions(k(ChosenOne),:) = LatticeCoords(ParticleLocation(k(ChosenOne)),:) + (-L/2 + (L/2 + L/2)*rand(1,2));
    else       
        for i = 0:(MDL-1)
            for j = (MDL+SDBase):(-1):(1+SDBase)
                newj = j + i*DL;
                newj1 = j-1 + i*DL;
                if j>DL
                    newj = j-DL + i*DL;
                end
                if (j-1) > DL
                    newj1 = j-1-DL + i*DL;
                end
                if N(newj1) ~= 0
                    k = find(ParticleLocation == newj1,N(newj1));
                    for p = 1:length(k)
                        ParticleLocation(k(p)) = newj;
                        N(newj1) = N(newj1) - 1;
                        N(newj) = N(newj) + 1;
                        Positions(k(p),:) = LatticeCoords(newj,:) + (-L/2 + (L/2 + L/2)*rand(1,2));
                    end
                end
%                 if j ~= (MDL+SDBase)
%                     N(i*DL + newj) = N(i*DL + newj1);
%                 else
%                     N(i*DL + newj) = N(i*DL + newj) + N(i*DL + newj1);
%                 end
            end
%             N(i*DL + SDBase) = 0;
        end
        if SDBase <= DL
            SDBase = SDBase + 1;
        else
            SDBase = SDBase + 1 - DL;
        end
    end
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

