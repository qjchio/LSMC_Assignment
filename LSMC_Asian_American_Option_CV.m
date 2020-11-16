function ratio = LSMC_Asian_American_Option_CV(K, S0, r, T, N, M)

%K is strike price
%S0 is intial price
%r is risk free rate
%T is Time to mat in years
%N is how many times you want to split T into
%M is number of sample paths

dt = T/N;
Y = zeros(M,N); 
S = S0 * ones(M, N+1);
dWt = sqrt(dt)*randn(M,N);
P = zeros(M,1);
%disp(dWt);

for i = 2:N+1
    S(:,i) = S(:,i-1) + ...
              r * S(:,i-1) * dt + ...
                 S(:,i-1) .* dWt(:,i-1) .* getVolatility(S(:,i-1),S0);   
end

S = S(:,2:N+1);
E = round(N/2);
A = zeros(M,N);
%Generate matrix A
for i = 1:M
    for j = E+1:N
        A(i,j) = 1/60 * sum(S(i,j-60: j-1));
    end
end
%disp(mean(A(:,N)));
for i = 1:M
   Y(i,N) = max(0, (K-A(i,N))); 
end
discount = exp(-r*dt);
for i = 1:M
    P(i,1) = max( K - S(i,N) ,0) * discount^N;
end
%backward Propagation
for k = N:-1:(E+2) %stop at 367, earliest we go in the loop is 366
    j = 0;
    for i = 1:M
        if(A(i,k-1) < K) %ITM 
            j = j+1;
            X2(j) = A(i,k-1);
            Y2(j) = discount*Y(i,k);
        end
    end
    p = polyfit(X2, Y2, 2);
    for i = 1:M
        if K - A(i, k-1) > polyval(p, A(i,k-1)) %exercise early
            Y(i,k-1) = max(K-A(i,k-1) ,0);
            %exEarly = exEarly+1;
        else
            Y(i,k-1) = discount * Y(i,k); %discount it backwards
        end
    end
end

W = exp(-r*dt*(E+1)) * Y(:,E+1);
U = P;
% disp(size(U));
% disp(size(W));
U_bar = mean(U);
cov_WU = cov(W,U);
c_opt = -cov_WU(1,1)/ var(U);
W_CV = W + c_opt * (U - U_bar);

meanPriceCV = mean(W_CV)
varPriceCV = var(W_CV)

% disp(size(W_CV))
% disp(size(U))

meanPriceCV1 = exp(-r*dt*(E+1)) * mean(Y(:,E+1))
varPriceCV1 = var(exp(-r*dt*(E+1)) *(Y(:,E+1)))

% U = exp(-r*dt*(E+1)) * (Y(:,E+1));
% W_CV(1:50,1)
% U(1:50,1)
% mean(W_CV)
% mean(U)

ratio = varPriceCV1/varPriceCV;

%
%plot price path
%
%time = 0:dt:T;
%disp(size(time));
%disp(size([S0 , S(1,:)]));
% for i = 1:1
%     plot(time, [S0, S(i,:)], 'b');
%     hold on
%     plot(time, [S0, S1(i,:)], 'r');
%     
% end
% yline(K);
% hold off
% title("Plot of Sample Paths Against Time")
% xlabel("Time in years")
% ylabel("Price")
% legend({'Stock Sample Price Path 1','Stock Sample Price Path 2'},'Location','northwest')
%
