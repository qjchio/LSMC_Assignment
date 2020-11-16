function [meanPrice, varPrice] = LSMC_Asian_Euro_Option(K, S0, r, T, N, M)
%K is strike price
%S0 is intial price
%r is risk free rate
%T is Time to mat in years
%N is how many times you want to split T into
%M is number of sample paths

dt = T/N;
Y = zeros(M,N); %N = 730 for now
S = S0 * ones(M, N+1);
%column runs from 1 to 730+1, first column is S0
dWt = sqrt(dt)*randn(M,N);
%disp(dWt);


for i = 2:N+1
    S(:,i) = S(:,i-1) + ...
              r * S(:,i-1) * dt + ...
                 S(:,i-1) .* dWt(:,i-1) .* getVolatility(S(:,i-1),S0);    
end

S = S(:,2:N+1);
%get rid of S0, so now M rows, N columns
%option path calculated
E = round(N/2);
A = zeros(M,N);
%M rows sample paths
%Number of columns is 365, 1st column is day 366 
%disp(mean(S(:,N)));
for i = 1:M
    for j = E+1:N
        %j goes from 1 to 365
        A(i,j) = 1/60 * sum(S(i,j-60: j-1));
        %disp(j-60);
        %disp(j-1);
    end
end

%disp(mean(A(:,N)));

for i = 1:M
   Y(i,N) = max(0, (K-A(i,N))); 
end

%disp(Y(:,N));

discount = exp(-r*dt);
meanPrice = mean(Y(:,N) * discount^N);
varPrice = var(Y(:,N) * discount^N);


[meanPrice, varPrice]


%
%plot price path
%
%time = 0:dt:T;
%disp(size(time));
%disp(size([S0 , S(1,:)]));
% for i = 1:M
%     plot(time, [S0, S(i,:)], 'b');
%     hold on
%     plot(time, [0 , A(i,:)], 'r');
% end
% yline(108);
% hold off
% title("Plot of Sample Paths Against Time")
% xlabel("Time in years")
% ylabel("Price")
% legend({'Stock Sample Price Path','Asian American Sample price Path'},'Location','northwest')

% 
% mean(S(:,N)*discount^N)
% mean(A(:,N)*discount^N)
% var(S(:,N)*discount^N)
% var(A(:,N)*discount^N)
