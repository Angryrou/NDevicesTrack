clear;
clc;
% Set total number of devices
N = 5;
% set sampling frequency
dt = 0.01;
% set total time of motion
T = 2 * pi;
% Time step set
t = 0 : dt : T;
% N is the total number of points
n = round(T/dt+1);
% set real trajectory
realP = zeros(2*N, n);
realV = zeros(2*N, n);
realA = zeros(2*N, n);
realD = 1;
for i = 1:N
    realP(2*i-1,:) = i * cos(t - pi/2);
    realP(2*i,:) = i * sin(t - pi/2);
    realV(2*i-1,:) = i * -sin(t - pi/2);
    realV(2*i,:) = i * cos(t - pi/2);
    realA(2*i-1,:) = i * -cos(t - pi/2);
    realA(2*i,:) = i * -sin(t - pi/2);
end
% set the process noise(m/s^2)
w = 0.4;
% set measure noise(m/s^2)
z = 0.1;
% set measured simulation acceleration
% rng('default');
measuredA = realA + (z + w) .* randn(size(realA));

%--  trajectory reconstruction with normal kalman filter  --%

% initialize state vector, distance estimation
stateofKF = zeros(6*N,n);
%dofKF = zeros(1,n);
% initiate state transition matrix (state = (ax1, vx1, px1, ay1, vy1, py1, ax2, vx2, px2, ay2, vy2, py2, ax3, vx3, px3, ay3, vy3, py3))
A = zeros(6*N, 6*N);
% initiate measurement matrix
C = zeros(2*N, 6*N);
% initiate process noise covariance
Q = zeros(6*N, 6*N);
% initiate measure noise covariance
R = zeros(2*N, 2*N);
% initiate forcast error covariance
P = zeros(6*N, 6*N);
for i = 1:2*N
    % set state transition matrix
    A(3*i-2:3*i, 3*i-2:3*i) = [1,0,0; dt,1,0; 1/2*dt*dt,dt,1];
    % set measurement matrix
    C(i, 3*i-2) = 1;
    % set process noise covariance
    Q(3*i-2, 3*i-2) = w^2;
    % set measure noise covariance
    R(i,i) = z^2;
    % set forcast error covariance
    P(3*i-2:3*i, 3*i-2:3*i) = 5 * diag([1,dt,dt*dt]);
end
% set system input
B = 0;
u = zeros(6*N, 1);

% do kalman filter
for i = 1 : n
    if i == 1
        for j = 1:N
            stateofKF(6*j-5,i) = measuredA(2*j-1,i);
            stateofKF(6*j-4,i) = j;
            stateofKF(6*j-3,i) = 0;
            stateofKF(6*j-2,i) = measuredA(2*j,i);
            stateofKF(6*j-1,i) = 0;
            stateofKF(6*j,i) = -j;
        end
    else
        [stateofKF(:,i),P] = kalmanFilter(stateofKF(:,i-1), P, A, B, u, C, measuredA(:,i), Q, R); 
    end
 %   dofKF(i) = sqrt((stateofKF(3,i) - stateofKF(9,i))^2 + (stateofKF(6,i) - stateofKF(12,i))^2);
end

pRMSofKF = zeros(2*N);
for i = 1:2*N
    pRMSofKF(i,:) = sqrt(sum((stateofKF(3*i,:) - realP(i,:)) .* (stateofKF(3*i,:) - realP(i,:))) / n);
end
%dRMSofKF = sqrt(sum((dofKF - realD) .* (dofKF - realD)) / N);

%-- trajectory reconstruction with equality-constrained MAUKF --%

stateofMAUKF = zeros(6*N,n);
%dofMAUKF = zeros(1,N);
% initiate process noise covariance
Q = zeros(6*N, 6*N);
% initiate measure noise covariance
R = zeros(3*N-1, 3*N-1);
% initiate forcast error covariance
P = zeros(6*N, 6*N);
for i = 1:2*N
    % set process noise covariance
    Q(3*i-2, 3*i-2) = w^2;
    % set measure noise covariance
    R(i,i) = z^2;
    % set forcast error covariance
    P(3*i-2:3*i, 3*i-2:3*i) = 5 * diag([1,dt,dt*dt]);
end
% set system input
u = zeros(6*N, 1);
% initiate MAUKF parameter
alpha = 1;
beta = 2;
kappa = 0;

for i = 1 : n
    if i == 1
        for j = 1:N
            stateofMAUKF(6*j-5,i) = measuredA(2*j-1,i);
            stateofMAUKF(6*j-4,i) = j;
            stateofMAUKF(6*j-3,i) = 0;
            stateofMAUKF(6*j-2,i) = measuredA(2*j,i);
            stateofMAUKF(6*j-1,i) = 0;
            stateofMAUKF(6*j,i) = -j;
        end
    else
        augmentedMeasure = zeros(3*N-1,1);
        augmentedMeasure(1:2*N,1) = measuredA(:,i);
        for k = 2*N+1 : 3*N-1
            augmentedMeasure(k,1) = realD;
        end
        [stateofMAUKF(:,i),P] = MAUKF(stateofMAUKF(:,i-1), P, u, dt, @ffun, @hfun, augmentedMeasure, Q, R, alpha, beta, kappa);
    end
%     dofMAUKF(i) = sqrt((stateofMAUKF(3,i) - stateofMAUKF(9,i))^2 + (stateofMAUKF(6,i) - stateofMAUKF(12,i))^2);
end
pRMSofMAUKF = zeros(2*N);
for i = 1:2*N
    pRMSofMAUKF(i,:) = sqrt(sum((stateofMAUKF(3*i,:) - realP(i,:)) .* (stateofMAUKF(3*i,:) - realP(i,:))) / n);
end
% dRMSofMAUKF = sqrt(sum((dofMAUKF - realD) .* (dofMAUKF - realD)) / N);

%-- consequence exhibitoin part --%
figure
hold on;
for i =1:N
    plot(realP(2*i-1,:),realP(2*i,:),'k');
    plot(stateofKF(6*i-3,:),stateofKF(6*i,:),'b');
    plot(stateofMAUKF(6*i-3,:),stateofMAUKF(6*i,:),'r');
end
hold off;
xlabel('x/m');
ylabel('y/m');
axis('square');
grid on;

