% (First) Zero Finding Method
% Should return the first zero of the  function collision(t+eps(t), tn(n), xn(:, n), v(:, n), k, epsilon, omega)

trigger = false;    % A boolean variable to terminate a while loop
safe = 0;           % A safety parameter to prevent an infinite while loop
safemax = 10^5;     % Max value for safe before the algorithm moves on.

% Initial conditons for time
t = tn(n);              % Time variable
tmax = t + max/mag(n);  % Maximum time until the ball is outside

% Inform me if the particle is outside the billaird
%  if xn(1, n)^2+xn(2, n)^2 > (1 + epsilon)^2
%      sprintf('THE BALL IS WAYYYY OUTSIDE')
%  elseif xn(1, n)^2+xn(2, n)^2 > (1 + epsilon*sin(omega*t)  )^2
%      sprintf('THE BALL IS  OUTSIDE')
%  end


% The main algorithm: 
% 1) Find a zero 
% 2) Find a minimun between the initial time t, and the new zero Z
% 3) If the minimun is less than zero, find a new zero (Z)
% 4) Repeat, until it found the first zero

Z = fzero(@(t)collision(t, tn(n), xn(:, n), v(:, n), k, epsilon, omega), [t + eps(t), tmax]);
while (trigger == false )
    
    min = fminbnd(@(a)collision(a, tn(n), xn(:, n), v(:, n), k, epsilon, omega), t + eps(t), Z);
   
    if collision(min, tn(n), xn(:, n), v(:, n), k, epsilon, omega) < 0 && t + eps(t) < min
        Z = fzero(@(t)collision(t, tn(n), xn(:, n), v(:, n), k, epsilon, omega), [t + eps(t), min]);
      else
        trigger = true;
    end
end

% Assign the time to that first zero
t = Z;

% Calculate position, et al...
xn(:, n+1) = xn(:, n) + v(:, n)*(t-tn(n));
    
% Calculate theta
if xn(1, n+1) > 0
    thetan(n+1)= atan(xn(2, n+1)/xn(1, n+1));
elseif xn(1, n+1) < 0
    thetan(n+1) = atan(xn(2, n+1)/xn(1, n+1)) + pi;
elseif ((xn(1, n+1)==0) && ( xn(2, n+1) > 0))
    thetan(n+1) = pi/2;
elseif ((xn(1, n+1) == 0) && ( xn(2, n+1) < 0))
    thetan(n+1) = 3*pi/2;
end
thetan(n+1) = mod(thetan(n+1), 2*pi);

% Calculate T, N, V
T(:, 1)= [ (R_0 ^ 2 + 0.2e1 * R_0 * epsilon * sin(omega * t) * sin(k * thetan(n+1)) + epsilon ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 - epsilon ^ 2 * cos(omega * t) ^ 2 + epsilon ^ 2 * cos(omega * t) ^ 2 * cos(k * thetan(n+1)) ^ 2 + epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 * cos(omega * t) ^ 2) ^ (-0.1e1 / 0.2e1) * (epsilon * sin(omega *  t) * cos(k * thetan(n+1)) * k * cos(thetan(n+1)) - sin(thetan(n+1)) * R_0 - sin(thetan(n+1)) * epsilon * sin(omega *  t) * sin(k * thetan(n+1))),  (R_0 ^ 2 + 0.2e1 * R_0 * epsilon * sin(omega * t) * sin(k * thetan(n+1)) + epsilon ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 - epsilon ^ 2 * cos(omega * t) ^ 2 + epsilon ^ 2 * cos(omega * t) ^ 2 * cos(k * thetan(n+1)) ^ 2 + epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 * cos(omega * t) ^ 2) ^ (-0.1e1 / 0.2e1) * (epsilon * sin(omega *  t) * cos(k * thetan(n+1)) * k * sin(thetan(n+1)) + cos(thetan(n+1)) * R_0 + cos(thetan(n+1)) * epsilon * sin(omega * t) * sin(k * thetan(n+1)))];
N(:, 1)= [ -T(2), T(1) ];
V(:, 1)= [epsilon * cos(omega * t) * omega * sin(k * thetan(n+1)) * cos(thetan(n+1)), epsilon * cos(omega * t) * omega * sin(k * thetan(n+1)) * sin(thetan(n+1))];

% Calcualte velocity
v(:, n+1)= v(:, n)-2*((N'*minus(v(:, n), V(:, 1)))*N);

% This loop corrects the solution found by fzero, and makes sure the
% solution stay inside the billaird 
% This is costly bevcause it repeatedly caclulates position and velocity
safe = 0;
while collision(t+eps(t), t+eps(t), xn(:, n+1), v(:, n+1), k, epsilon, omega) <= 0 && safe <= safemax
    
    % Take a time step backwards, and recalculate all essential variables
    t = t - eps(t);
    xn(:, n+1) = xn(:, n)+v(:, n)*(t-tn(n));
    
    % Calculate theta
    if xn(1, n+1) > 0
        thetan(n+1) = atan(xn(2, n+1)/xn(1, n+1));
    elseif xn(1, n+1) < 0
        thetan(n+1) = atan(xn(2, n+1)/xn(1, n+1)) + pi;
    elseif ((xn(1, n+1) == 0) && ( xn(2, n+1) > 0))
        thetan(n+1) = pi/2;
    elseif ((xn(1, n+1) == 0) && ( xn(2, n+1) < 0))
        thetan(n+1) = 3*pi/2;
    end
    thetan(n+1) = mod(thetan(n+1), 2*pi);

    % Calculate T, N, V
    T(:, 1) = [ (R_0 ^ 2 + 0.2e1 * R_0 * epsilon * sin(omega * t) * sin(k * thetan(n+1)) + epsilon ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 - epsilon ^ 2 * cos(omega * t) ^ 2 + epsilon ^ 2 * cos(omega * t) ^ 2 * cos(k * thetan(n+1)) ^ 2 + epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 * cos(omega * t) ^ 2) ^ (-0.1e1 / 0.2e1) * (epsilon * sin(omega *  t) * cos(k * thetan(n+1)) * k * cos(thetan(n+1)) - sin(thetan(n+1)) * R_0 - sin(thetan(n+1)) * epsilon * sin(omega *  t) * sin(k * thetan(n+1))),  (R_0 ^ 2 + 0.2e1 * R_0 * epsilon * sin(omega * t) * sin(k * thetan(n+1)) + epsilon ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 - epsilon ^ 2 * cos(omega * t) ^ 2 + epsilon ^ 2 * cos(omega * t) ^ 2 * cos(k * thetan(n+1)) ^ 2 + epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 - epsilon ^ 2 * cos(k * thetan(n+1)) ^ 2 * k ^ 2 * cos(omega * t) ^ 2) ^ (-0.1e1 / 0.2e1) * (epsilon * sin(omega *  t) * cos(k * thetan(n+1)) * k * sin(thetan(n+1)) + cos(thetan(n+1)) * R_0 + cos(thetan(n+1)) * epsilon * sin(omega * t) * sin(k * thetan(n+1)))];
    N(:, 1) = [ -T(2), T(1) ];
    V(:, 1) = [epsilon * cos(omega * t) * omega * sin(k * thetan(n+1)) * cos(thetan(n+1)), epsilon * cos(omega * t) * omega * sin(k * thetan(n+1)) * sin(thetan(n+1))];

    % Calcualte velocity
    v(:, n+1) = v(:, n)-2*( (N'*minus(v(:, n), V(:, 1)))*N);
 
    if safe > 1
        safe
    end
    safe = safe + 1;
end

% Assign the new collision time
tn(n+1) = t;
    
% Calculate the magnitude of the velocity
 mag(1, n+1) = sqrt(v(:, n+1)'* v(:, n+1));
    
% Calculate alpha
if (v(:, n+1)'*N >0)
    alpha(n+1) = acos((v(:, n+1)'*T)/mag(1, n+1));
else
    alpha(n+1) = 2*pi - acos((v(:, n+1)'*T)/mag(1, n+1));
end

% Calculate distance bewtwen last two collisions
dist(n+1) = xn(:, n)'*xn(:, n+1);
    
% Calculate Modded Time 
modt(n+1) = mod(tn(n+1), 2*pi);
    
% Maxv
if v(n+1) > vmax
    vmax = v(n+1);
end

