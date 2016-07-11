function [newX,newP] = MAUKF(oldX,oldP,u,dt,ffun,hfun,newY,Q,R,alpha,beta,kappa)
    % calculate sigma point matrix and weight
    [xPts,wPts,nPts] = sigmaPoints(oldX,oldP,alpha,beta,kappa);
    sigmaM = wPts(1:nPts);
    sigmaC = wPts(1:nPts);
    sigmaC(1) = wPts(nPts + 1);
    xPtsEst = zeros(size(xPts));
    for i = 1 : nPts
        xPtsEst(:,i) = ffun(xPts(:,i),u,dt);
    end
    xEst = zeros(size(oldX));
    PxxEst = zeros(size(oldP));
    for i = 1 : nPts
        xEst = xEst + sigmaM(i) .* xPtsEst(:,i);
    end
    for i = 1 : nPts
        PxxEst = PxxEst + sigmaC(i) .* (xPtsEst(:,i) - xEst) * (xPtsEst(:,i) - xEst)';
    end
    PxxEst = PxxEst + Q;
    
    yPtsEst = zeros(size(newY,1),nPts);
    for i = 1 : nPts
        yPtsEst(:,i) = hfun(xPtsEst(:,i));
    end
    yEst = zeros(size(newY));
    PyyEst = zeros(size(newY,1),size(newY,1));
    PxyEst = zeros(size(oldX,1),size(newY,1));
    for i = 1 : nPts
        yEst = yEst + sigmaM(i) .* yPtsEst(:,i);
    end
    for i = 1 : nPts
        PyyEst = PyyEst + sigmaC(i) .* (yPtsEst(:,i) - yEst) * (yPtsEst(:,i) - yEst)';
    end
    PyyEst = PyyEst + R;
    for i = 1 : nPts
        PxyEst = PxyEst + sigmaC(i) .* (xPtsEst(:,i) - xEst) * (yPtsEst(:,i) - yEst)';
    end
    Kg = PxyEst / PyyEst;
    newX = xEst + Kg * (newY - yEst);
    newP = PxxEst - Kg * PyyEst * Kg';
end