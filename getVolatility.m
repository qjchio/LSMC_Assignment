function sigma = getVolatility(X,S)

first = X <= 100;
second = X > 100;

bool1 = first;
bool2 = second;

first = first .* X;
first = 0.25 + 0.02 * ( 1 - first/S);
first = first .* bool1;

second = second .* X;
second = max(0.001, 0.25 - 0.01 * ( second/S - 1));
second = second .* bool2;

sigma = first + second;


end