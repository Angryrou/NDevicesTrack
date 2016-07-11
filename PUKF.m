function [newX,newP] = PUKF(oldX,oldP,u,dt,ffun,hfun,gfun,newY,d,Q,R,alpha,beta,kappa)
    % calculate sigma point matrix and weight
    [xPts,wPts,nPts] = sigmaPoints(oldX,oldP,alpha,beta,kappa);
    sigmaM = wPts(1:nPts);
    sigmaC = wPts(1:nPts);
    sigmaC(1) = wPts(nPts + 1);
    % project state vector according to equality-constraint
    dPts = zeros(1,nPts);
    for i = 1 : nPts
        dPts(i) = gfun(xPts(:,i));
    end
    dEst = sum(sigmaM .* dPts);
    Pdd = 0;
    for i = 1 : nPts
        Pdd = Pdd + sigmaC(i) * (dPts(i) - dEst)^2;
    end
    Pxd = zeros(size(oldX,1),1);
    for i = 1 : nPts
        Pxd = Pxd + sigmaC(i) .* (xPts(:,i) - oldX) * (dPts(i) - dEst)';
    end
    Kgp = Pxd / Pdd;
    projectedX = oldX + Kgp * (d - dEst);
    Pxxp = oldP - Kgp * Pdd * Kgp';
    
    [xPts,wPts,nPts] = sigmaPoints(projectedX,Pxxp,alpha,beta,kappa);
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