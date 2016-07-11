function yEst = hfun(x)
    N = size(x,1)/6;
    % initiate measurement matrix
    C = zeros(2*N, 6*N);
    for i = 1:2*N
        % set measurement matrix
        C(i, 3*i-2) = 1;
    end
    d = zeros(N-1,1);
    for i = 1:N-1
        d(i) = sqrt((x(6*i+3) - x(6*i-3))^2 + (x(6*i+6) - x(6*i))^2);
    end
    yEst = [C * x;d];
end