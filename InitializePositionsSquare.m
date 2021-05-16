function Positions = InitializePositionsSquare(K,L)
N =  round(sqrt(K));
Positions = zeros(K,2);
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Positions(1:N,:) = Positions(1:N,:) + [0,N-0.5];
Positions(1:N,1) = (1:N)' - 0.5;

for i = 2:N
    Positions((1:N) + (i-1)*N,1) = Positions(1:N,1);
    Positions((1:N) + (i-1)*N,2) = Positions(1:N,2) - (i-1);
end

Positions = Positions + 0*[-N/2,-N/2];
Positions = Positions*L;

end


