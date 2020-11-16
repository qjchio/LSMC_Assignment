function ratio = LSMC_Asian_American_Option_IS(K, S0, r, T, N, M)

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
Y1 = zeros(M,N);
S1 = S0 * ones(M, N+1);
%disp(dWt);
%Change r_1 here v
r1 = r - 0.05;
RN = zeros(M,N);

for i = 2:N+1
    vol = getVolatility(S(:,i-1),S0);
    S(:,i) = S(:,i-1) + ...
              r * S(:,i-1) * dt + ...
                 S(:,i-1) .* dWt(:,i-1) .* vol;  
     
    vol1 = getVolatility(S1(:,i-1),S0);
    
    theta = (r-r1)/vol;
    
    S1(:,i) = S1(:,i-1) + ...
              r1 * S1(:,i-1) * dt + ...
                 S1(:,i-1) .* dWt(:,i-1) .* vol1;
%     disp(size(theta));
%     disp(size(dWt(:,i-1)));
%     disp(size(vol1));
%     disp(size(RN(:,i-1)));
%     disp(mean(S(:,i)-S1(:,i)));
      
    theta = theta';
   RN(:,i-1) = exp( - theta.^2 * dt / 2 + theta .* dWt(:,i-1));
end

S = S(:,2:N+1);
S1 = S1(:,2:N+1);
E = round(N/2);
A = zeros(M,N);
A1 = zeros(M,N);

for i = 1:M
    for j = E+1:N
        A(i,j) = 1/60 * sum(S(i,j-60: j-1));
        A1(i,j) = 1/60 * sum(S1(i,j-60:j-1));
    end
end

%disp(mean(A(:,N)));

for i = 1:M
   Y(i,N) = max(0, (K-A(i,N))); 
   Y1(i,N) = max(0, (K-A1(i,N)));
end
% sum((Y(:,N) > 0))
% sum((Y1(:,N) > 0))

discount = exp(-r*dt);

for k = N:-1:(E+2) %stop at 367, earliest we go in the loop is 366
    j = 0;
    j1 = 0;
    for i = 1:M
        if(A(i,k-1) < K) %ITM 
            j = j+1;
            X2(j) = A(i,k-1);
            Y2(j) = discount*Y(i,k);
        end
        if(A1(i,k-1) < K) %ITM 
            j1 = j1+1;
            X3(j1) = A1(i,k-1);
            Y3(j1) = discount*Y1(i,k);
        end
    end
    
    p = polyfit(X2, Y2, 2);
    p1 = polyfit(X3, Y3, 2);
    for i = 1:M
        if K - A(i, k-1) > polyval(p, A(i,k-1)) %exercise early
            Y(i,k-1) = max(K-A(i,k-1) ,0);
            
        else
            Y(i,k-1) = discount * Y(i,k); %discount it backwards
        end
        if K - A1(i, k-1) > polyval(p1, A1(i,k-1)) %exercise early
            Y1(i,k-1) = max(K-A1(i,k-1) ,0);
            
        else
            Y1(i,k-1) = discount * Y1(i,k); %discount it backwards
        end
    end
end


%disp(mean(Y(:,N)));
%disp(exp(-r*dt*(E+1)));

RN = prod(RN,2);

out = exp(-r*dt*(E+1)) * Y(:,E+1);

meanPriceIS = mean(out)
varPriceIS = var(out)

out1 = exp(-r*dt*(E+1)) * Y1(:,E+1) .* RN;
meanPriceIS1 = mean(out1)
varPriceIS1 = var(out1)

%disp(sdPrice);
% sum(out > 0)
% sum(out1 > 0)


ratio = varPriceIS/varPriceIS1;




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


% mean(S(:,N)*discount^N)
% mean(A(:,N)*discount^N)
% var(S(:,N)*discount^N)
% var(A(:,N)*discount^N)
