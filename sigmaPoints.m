function [xPts,wPts,nPts] = sigmaPoints(x,P,alpha,beta,kappa)
    n = size(x,1);
    nPts = 2 * n +1;
    %calculate lamda according to scaling parameters
    lamda = alpha^2 * (n + kappa) - n;
    %calculate sigma point matrix and the weight on the points
    xPts = [zeros(size(P,1),1), -(chol((n + lamda) * P))', (chol((n + lamda) * P))'];
    xPts = xPts + repmat(x,1,nPts);
    wPts = [lamda, 0.5 * ones(1,nPts-1), 0]/(n + lamda);
    wPts(nPts+1) = wPts(1) + 1 - alpha^2 + beta;
end