function xEst = ffun(x,u,dt)
    N = size(x,1)/6;
    A = zeros(6*N, 6*N);
    for i = 1:2*N
        % set state transition matrix
        A(3*i-2:3*i, 3*i-2:3*i) = [1,0,0; dt,1,0; 1/2*dt*dt,dt,1];
    end
    xEst = A * x + u;
end
    