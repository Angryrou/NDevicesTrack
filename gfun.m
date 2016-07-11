function d = gfun(x)
    d = sqrt((x(3) - x(9))^2 + (x(6) - x(12))^2);
    % or x' * D * x where D = [diag([0,0,1,0,0,1]),diag([0,0,-1,0,0,-1]);diag([0,0,-1,0,0,-1]),diag([0,0,1,0,0,1])]
end

