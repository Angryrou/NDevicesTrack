function [newX,newP] = kalmanFilter(oldX, oldP, A, B, u, C, newY, Q, R)
    xEst = A * oldX + B .* u;
    pEst = A * oldP * A' + Q;
    Kg = pEst * C' / (C * pEst * C' + R);
    newX = xEst + Kg * (newY - C * xEst);
    [m , ~] = size(pEst);
    newP = (eye(m) - Kg * C) * pEst;
end 