function Positions = InitializePositionsCube(K,L)
N =  round(K^(1/3));
Positions = zeros(K,3);
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Positions(1:N,:) = Positions(1:N,:) + [0,N-0.5,0.5];
Positions(1:N,1) = (1:N)' - 0.5;

for i = 2:N
    Positions((1:N) + (i-1)*N,[1,3]) = Positions(1:N,[1,3]);
    Positions((1:N) + (i-1)*N,2) = Positions(1:N,2) - (i-1);
end

for k = 2:N
    Positions((1:(N*N)) + (k-1)*N*N,1:2) = Positions(1:(N*N),1:2);
    Positions((1:(N*N)) + (k-1)*N*N,3) = Positions(1:(N*N),3) + (k-1);
end

Positions = Positions + [-N/2,-N/2,-N/2];
Positions = Positions*L;

end


