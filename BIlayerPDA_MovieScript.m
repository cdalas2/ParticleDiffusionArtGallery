clear all
close all
clc
DIM = 2;
HeadOutlinePositionsData = load("BilayerPDAHeadOutline.o");
HeadOutlineSender = HeadOutlinePositionsData(1:end-1,1)+1;
HeadOutlineReceiver = HeadOutlinePositionsData(1:end-1,2)+1;
HeadOutlinerho = HeadOutlinePositionsData(1:end-1,3);
HeadOutlinesimT = HeadOutlinePositionsData(1:end-1,4);
HeadOutlineK = HeadOutlinePositionsData(end,1);

% SquareN = ones(SquareK,1);
HeadOutlineN = zeros(HeadOutlineK,1); HeadOutlineN(HeadOutlineSender(1)) = HeadOutlineK;
HeadOutlinerho_max = HeadOutlinerho(end)
HeadOutlineLCELLS_PER_LENGTH_SCALE = HeadOutlinePositionsData(end,2);
HeadOutlineTIME_MAX = HeadOutlinePositionsData(end,4);
HeadOutlineLENGTH_SCALE = 10; %1 micrometer
HeadOutlineL = HeadOutlineLENGTH_SCALE/HeadOutlineLCELLS_PER_LENGTH_SCALE;
HeadOutlineLatticeCoords = InitializePositionsSquare(1600,10);
HeadOutlinePositions = zeros(HeadOutlineN(HeadOutlineSender(1)),DIM);
for i = 1:HeadOutlineN(HeadOutlineSender(1))
    %     SquarePositions(i,:)= SquareLatticeCoords(i,:) + (-SquareL/2 + (SquareL/2 + SquareL/2)*rand(1,2));
    HeadOutlinePositions(i,:)= HeadOutlineLatticeCoords(HeadOutlineSender(1),:) + (-HeadOutlineL/2 + (HeadOutlineL/2 + HeadOutlineL/2)*rand(1,2));
end
% SquareParticleLocation = 1:SquareK;
HeadOutlineParticleLocation = ones(HeadOutlineN(HeadOutlineSender(1)),1)*HeadOutlineSender(1);


% % % Din = 20;
% % % Dout = 1000;

HeadInsidePositionsData = load("BilayerPDAHeadInside.o");
HeadInsideSender = HeadInsidePositionsData(1:end-1,1)+1;
HeadInsideReceiver = HeadInsidePositionsData(1:end-1,2)+1;
HeadInsiderho = HeadInsidePositionsData(1:end-1,3);
HeadInsidesimT = HeadInsidePositionsData(1:end-1,4);
HeadInsideK = HeadInsidePositionsData(end,1);
% N = ones(K,1);
HeadInsideN = zeros(HeadInsideK,1); HeadInsideN(HeadInsideSender(1)) = HeadInsideK;
HeadInsiderho_max = HeadInsiderho(end)
HeadInsideLCELLS_PER_LENGTH_SCALE = HeadInsidePositionsData(end,2);
HeadInsideTIME_MAX = HeadInsidePositionsData(end,4);
HeadInsideLENGTH_SCALE = 10; %1 micrometer
HeadInsideL = HeadInsideLENGTH_SCALE/HeadInsideLCELLS_PER_LENGTH_SCALE;
HeadInsideLatticeCoords = InitializePositionsSquare(1600,10);
HeadInsidePositions = zeros(HeadInsideN(HeadInsideSender(1)),DIM);
for i = 1:HeadInsideN(HeadInsideSender(1))
    %     Positions(i,:)= LatticeCoords(Sender(1),:) + (-L/2 + (L/2 + L/2)*rand(1,3));
    HeadInsidePositions(i,:)= HeadInsideLatticeCoords(HeadInsideSender(1),:) + (-HeadInsideL/2 + (HeadInsideL/2 + HeadInsideL/2)*rand(1,2));
    
end
% ParticleLocation = (1:K)';
% SDBase = 1;
DL = HeadInsideK^(1/2);
% MDL = 18;
HeadInsideParticleLocation = ones(HeadInsideN(HeadInsideSender(1)),1)*HeadInsideSender(1);


ChainsPositionsData = load("BilayerPDAChains.o");
ChainsSender = ChainsPositionsData(1:end-1,1)+1;
ChainsReceiver = ChainsPositionsData(1:end-1,2)+1;
Chainsrho = ChainsPositionsData(1:end-1,3);
ChainssimT = ChainsPositionsData(1:end-1,4);
ChainsK = ChainsPositionsData(end,1);

% ChainsN = ones(ChainsK,1);
ChainsN = zeros(ChainsK,1); ChainsN(ChainsSender(1)) = ChainsK;
Chainsrho_max = Chainsrho(end)
ChainsLCELLS_PER_LENGTH_SCALE = ChainsPositionsData(end,2);
ChainsTIME_MAX = ChainsPositionsData(end,4);
ChainsLENGTH_SCALE = 10; %1 micrometer
ChainsL = ChainsLENGTH_SCALE/ChainsLCELLS_PER_LENGTH_SCALE;
ChainsLatticeCoords = InitializePositionsSquare(1600,10);
ChainsPositions = zeros(ChainsN(ChainsSender(1)),DIM);
for i = 1:ChainsN(ChainsSender(1))
    %     ChainsPositions(i,:)= ChainsLatticeCoords(i,:) + (-ChainsL/2 + (ChainsL/2 + ChainsL/2)*rand(1,2));
    ChainsPositions(i,:)= ChainsLatticeCoords(ChainsSender(1),:) + (-ChainsL/2 + (ChainsL/2 + ChainsL/2)*rand(1,2));
end
% ChainsParticleLocation = 1:ChainsK;
ChainsParticleLocation = ones(ChainsN(ChainsSender(1)),1)*ChainsSender(1);


EmptySpaceInsidePositionsData = load("BilayerPDAEmptySpaceInside.o");
EmptySpaceInsideSender = EmptySpaceInsidePositionsData(1:end-1,1)+1;
EmptySpaceInsideReceiver = EmptySpaceInsidePositionsData(1:end-1,2)+1;
EmptySpaceInsiderho = EmptySpaceInsidePositionsData(1:end-1,3);
EmptySpaceInsidesimT = EmptySpaceInsidePositionsData(1:end-1,4);
EmptySpaceInsideK = EmptySpaceInsidePositionsData(end,1);

% EmptySpaceInsideN = ones(EmptySpaceInsideK,1);
EmptySpaceInsideN = zeros(EmptySpaceInsideK,1); EmptySpaceInsideN(EmptySpaceInsideSender(1)) = EmptySpaceInsideK;
EmptySpaceInsiderho_max = EmptySpaceInsiderho(end)
EmptySpaceInsideLCELLS_PER_LENGTH_SCALE = EmptySpaceInsidePositionsData(end,2);
EmptySpaceInsideTIME_MAX = EmptySpaceInsidePositionsData(end,4);
EmptySpaceInsideLENGTH_SCALE = 10; %1 micrometer
EmptySpaceInsideL = EmptySpaceInsideLENGTH_SCALE/EmptySpaceInsideLCELLS_PER_LENGTH_SCALE;
EmptySpaceInsideLatticeCoords = InitializePositionsSquare(1600,10);
EmptySpaceInsidePositions = zeros(EmptySpaceInsideN(EmptySpaceInsideSender(1)),DIM);
for i = 1:EmptySpaceInsideN(EmptySpaceInsideSender(1))
    %     EmptySpaceInsidePositions(i,:)= EmptySpaceInsideLatticeCoords(i,:) + (-EmptySpaceInsideL/2 + (EmptySpaceInsideL/2 + EmptySpaceInsideL/2)*rand(1,2));
    EmptySpaceInsidePositions(i,:)= EmptySpaceInsideLatticeCoords(EmptySpaceInsideSender(1),:) + (-EmptySpaceInsideL/2 + (EmptySpaceInsideL/2 + EmptySpaceInsideL/2)*rand(1,2));
end
% EmptySpaceInsideParticleLocation = 1:EmptySpaceInsideK;
EmptySpaceInsideParticleLocation = ones(EmptySpaceInsideN(EmptySpaceInsideSender(1)),1)*EmptySpaceInsideSender(1);


figure(34)
hold on
xlabel('$x$ ($\mu$m)','interpreter','latex','fontsize',19);
ylabel('$y$ ($\mu$m)','interpreter','latex','fontsize',19);
set(gca,'XColor', 'none','YColor','none')
aa = DL * HeadInsideL/2;
set(gca, 'XLim', [0 2*aa], 'YLim', [0 2*aa]);
% set(gca,'Color','none');
set(gca,'Color',[0.00,0.45,0.74]);
ms = 17;
% hp=plot(Positions(:,1),Positions(:,2),'.','MarkerSize',30,'MarkerEdgeColor',[255,204,0]/255,...
%     'MarkerFaceColor',[255,204,0]/255);
% HeadOutlinehp=plot(HeadOutlinePositions(:,1),HeadOutlinePositions(:,2),'.','MarkerSize',15,'MarkerEdgeColor',[0.85,0.33,0.10],...
%     'MarkerFaceColor',[0.85,0.33,0.10]);
HeadOutlinehp=plot(HeadOutlinePositions(:,1),HeadOutlinePositions(:,2),'.','MarkerSize',ms,'MarkerEdgeColor',[0.93,0.69,0.13],...
    'MarkerFaceColor',[0.93,0.69,0.13]);
HeadInsidehp=plot(HeadInsidePositions(:,1),HeadInsidePositions(:,2),'.','MarkerSize',ms,'MarkerEdgeColor',[0.93,0.69,0.13],...
    'MarkerFaceColor',[0.93,0.69,0.13]);
Chainshp=plot(ChainsPositions(:,1),ChainsPositions(:,2),'.','MarkerSize',ms,'MarkerEdgeColor',[0.60,0.20,0.48],...
    'MarkerFaceColor',[0.60,0.20,0.48]);
EmptySpaceInsidehp=plot(EmptySpaceInsidePositions(:,1),EmptySpaceInsidePositions(:,2),'.','MarkerSize',ms,'MarkerEdgeColor',[0.84,0.91,0.85],...
    'MarkerFaceColor',[0.84,0.91,0.85]);

% comphp=plot([SquarePositions(:,1);Positions(:,1)],[SquarePositions(:,2);Positions(:,2)],'r.','MarkerSize',30);
v = VideoWriter('BilayerPDA');
v.Quality = 95;
open(v);
M(1) = getframe(gcf);
writeVideo(v,M(1));
Tsteps = length(HeadInsidesimT);

for step = 1:Tsteps
    %     pause((simT(i+1)-simT(i))/100);
    HeadOutlineChosenOne = randi(HeadOutlineN(HeadOutlineSender(step)));
    HeadOutlinek = find(HeadOutlineParticleLocation == HeadOutlineSender(step),HeadOutlineN(HeadOutlineSender(step)));
    HeadOutlineParticleLocation(HeadOutlinek(HeadOutlineChosenOne)) = HeadOutlineReceiver(step);
    HeadOutlineN(HeadOutlineSender(step)) = HeadOutlineN(HeadOutlineSender(step)) - 1;
    HeadOutlineN(HeadOutlineReceiver(step)) = HeadOutlineN(HeadOutlineReceiver(step)) + 1;
    HeadOutlinePositions(HeadOutlinek(HeadOutlineChosenOne),:) = HeadOutlineLatticeCoords(HeadOutlineReceiver(step),:) + (-HeadOutlineL/2 + (HeadOutlineL/2 + HeadOutlineL/2)*rand(1,2));
    
    HeadInsideChosenOne = randi(HeadInsideN(HeadInsideSender(step)));
    HeadInsidek = find(HeadInsideParticleLocation == HeadInsideSender(step),HeadInsideN(HeadInsideSender(step)));
    HeadInsideParticleLocation(HeadInsidek(HeadInsideChosenOne)) = HeadInsideReceiver(step);
    HeadInsideN(HeadInsideSender(step)) = HeadInsideN(HeadInsideSender(step)) - 1;
    HeadInsideN(HeadInsideReceiver(step)) = HeadInsideN(HeadInsideReceiver(step)) + 1;
    HeadInsidePositions(HeadInsidek(HeadInsideChosenOne),:) = HeadInsideLatticeCoords(HeadInsideParticleLocation(HeadInsidek(HeadInsideChosenOne)),:) + (-HeadInsideL/2 + (HeadInsideL/2 + HeadInsideL/2)*rand(1,2));
    
    ChainsChosenOne = randi(ChainsN(ChainsSender(step)));
    Chainsk = find(ChainsParticleLocation == ChainsSender(step),ChainsN(ChainsSender(step)));
    ChainsParticleLocation(Chainsk(ChainsChosenOne)) = ChainsReceiver(step);
    ChainsN(ChainsSender(step)) = ChainsN(ChainsSender(step)) - 1;
    ChainsN(ChainsReceiver(step)) = ChainsN(ChainsReceiver(step)) + 1;
    ChainsPositions(Chainsk(ChainsChosenOne),:) = ChainsLatticeCoords(ChainsReceiver(step),:) + (-ChainsL/2 + (ChainsL/2 + ChainsL/2)*rand(1,2));
    
    EmptySpaceInsideChosenOne = randi(EmptySpaceInsideN(EmptySpaceInsideSender(step)));
    EmptySpaceInsidek = find(EmptySpaceInsideParticleLocation == EmptySpaceInsideSender(step),EmptySpaceInsideN(EmptySpaceInsideSender(step)));
    EmptySpaceInsideParticleLocation(EmptySpaceInsidek(EmptySpaceInsideChosenOne)) = EmptySpaceInsideReceiver(step);
    EmptySpaceInsideN(EmptySpaceInsideSender(step)) = EmptySpaceInsideN(EmptySpaceInsideSender(step)) - 1;
    EmptySpaceInsideN(EmptySpaceInsideReceiver(step)) = EmptySpaceInsideN(EmptySpaceInsideReceiver(step)) + 1;
    EmptySpaceInsidePositions(EmptySpaceInsidek(EmptySpaceInsideChosenOne),:) = EmptySpaceInsideLatticeCoords(EmptySpaceInsideReceiver(step),:) + (-EmptySpaceInsideL/2 + (EmptySpaceInsideL/2 + EmptySpaceInsideL/2)*rand(1,2));
    
    
    set(HeadInsidehp,'XData',HeadInsidePositions(:,1),'YData',HeadInsidePositions(:,2));
%     if mod(step,5000) == 0
%         drawnow;
%     end
    %     set(comphp,'XData',[SquarePositions(:,1);Positions(:,1)],'YData',[SquarePositions(:,2);Positions(:,2)]);
    set(HeadOutlinehp,'XData',HeadOutlinePositions(:,1),'YData',HeadOutlinePositions(:,2));
    set(Chainshp,'XData',ChainsPositions(:,1),'YData',ChainsPositions(:,2));
    set(EmptySpaceInsidehp,'XData',EmptySpaceInsidePositions(:,1),'YData',EmptySpaceInsidePositions(:,2));
    if mod(step,1200) == 0
        drawnow;
        HeadInsidesimT(step)
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

