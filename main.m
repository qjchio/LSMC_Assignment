tic
K = 108;
S0 = 100;
r = 0.08;
T = 2;
N = 365*2;%by constraint of the question, leave N as 365*2 because A_n is defined as average of past 60 days. //for simplicity
M = 10000;

%%Comment out respective lines to get result

%Euro version just to check if prices are about the same
%[meanPrice, varPrice] = LSMC_Asian_Euro_Option(K,S0,r,T,N,M)

%Asian American put option pricing
[meanPrice, varPrice] = LSMC_Asian_American_Option(K,S0,r,T,N,M)

%AntiThetic Variate Method
ratioAV = LSMC_Asian_American_Option_AV(K,S0,r,T,N,M)

%Control Variate Method
ratioCV = LSMC_Asian_American_Option_CV(K,S0,r,T,N,M)

%Importance Sampling Method
ratioIS = LSMC_Asian_American_Option_IS(K,S0,r,T,N,M)

toc

%about 130 seconds