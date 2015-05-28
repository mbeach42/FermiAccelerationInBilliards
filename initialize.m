% Initialize.m: This is where the entire project starts
% In this file I initialize initial parameters, arrays, and start the
% simulation of a ball udner going M collision with a moving boundary
% perturbed by 

tic     % How long does the simulation take?

% Turn off LaTeX interpreter for making plots
% set(0, 'defaulttextinterpreter', 'none');

% Fundamental parameters, 
R_0 = 1;                    % R_0 = radius of the circle, 
omega = 1;                  % omega = the angular frequency, 
k = 2;                      % k = the wavenumber, 
M = 100;                    % M = number of collisions to be simulated, M>50 for animation to work

% Vectors
nn = zeros(1, M);           % Store the current n value
v = zeros(2, M);            % Velocity of particle
V = zeros(2, 1);            % Boundry velocity vectors
T = zeros(2, 1);            % Boundry Tangent vector
N = zeros(2, 1);            % Boundry Normal vector
xn = zeros(2, M);           % Cartesian position of particle
mag = zeros(1, M);          % Magnitude of velocity
thetan = zeros(1, M);       % Theta, angular direction
tn = zeros(1, M);           % Collision times
alpha = zeros(1, M);        % Angle leaving the collision
dist = zeros(1, M);         % Distance between collision
modt = zeros(1, M);         % Collision time modded by 2*Pi (Phase)

% Initial conditions for time, modded time and the maximun velocity
tn(1) = 2*pi;               % The 0th collision has time 2*pi for convenience
modt(1) = 0; 
vmax = 0;
mag(1) = 1;                 % The initial velocity is always 1
nn(1) = 1;                  % The initial collision number is 1

epsilonrange = 9:9;         % The range of the perturbation amplitude (epsilon)   
rangelength = length(epsilonrange);     % The number of epsilons used

% Initializing the linear and exponential fitting vectors.
slope = zeros(rangelength, 5, 5);
avgslope = zeros(rangelength, 1);
linslope = zeros(rangelength, 5, 5);
avglinslope = zeros(rangelength, 1);
lastvarray = zeros(rangelength, 5, 5);
lastvavg = zeros(rangelength, 1);

% Close all figures
close all

% Initialize a counting variable
counter = 0;

% The initial conditions that need to be iterated over
% For k = 1, do epsilon = 0 .. 1, theta = -pi/2 .. pi/2, alpha = -pi/2 .. pi/2
% For k = 2, do epsilon = 0 .. 1, theta = 0 .. pi/2, alpha = 0 .. pi/2

% This is the main loop which iterates over initial conditions
% It does so in 3 loops:
    % -range over epsilon (epsilonrange)
    % -range over initial theta position of partilce (thetarange)
    % -range over initial leaving trajectory angle alpha (alpharange)
    
for h = epsilonrange
    
    epsilon = h/10 - 1/10;    % This actually changes the epsilon parameter
    max  = 2 + 2*epsilon;     % This is implimented in 'zerofindingmethod.m' as the maximum range
    
    % Usually thetarange = 1:5
    for thetarange = 2:2
        
        thetan(1) = (thetarange - 1)*pi/8;   % This actually changes the theta parameter
        
        % Calculate the position at theta, 
        xn(1, 1) = cos(thetan(1));
        xn(2, 1) = sin(thetan(1));

        % Usually alpharange = 2:5
        for alpharange = 2:2
            
            % To save the dataset
            % counter = counter + 1;
            % savefile = sprintf('DataSet_% d_% d_% d_% d.mat', k, epsilon, thetarange, alpharange);

            alpha(1) = (alpharange-1)*pi/8 ;  % This actually changes the alpha parameter
            
            % Calculate the initial velocity vector
            v(:, 1) = [cos(thetan(1) + alpha(1) + pi/2), sin(thetan(1) + alpha(1) + pi/2) ];
            
            % Notify me that it's beginning the simulation
            sprintf('Simulating...')

            % Iterate over M collisions
            for n = 1:(M-1)
                nn(n+1) = n+1;        % Record the collision number
                zerofindingmethod   % Call the method that finds the collision time, et al.
            end
            
            % Notify me that it's done simulating and is now calculating the
            sprintf('Calculating fits...')
            
            % Store the final velocity
            lastvarray(h, thetarange, alpharange) = mag(n+1);
            
            % Find a power law (n^a) fit for velocity vs time
            X = log(nn);
            Y = log(mag);
            p = polyfit(X, Y, 1);
            pval = polyval(p, X);
            slope(h, thetarange, alpharange) = p(1)*M/tn(M);
            
            % Find a linear fit for the velocity vs time
            ptwo = polyfit(nn, mag, 1);
            pvaltwo = polyval(ptwo, nn);
            linslope(h, thetarange, alpharange) = ptwo(1)*M/tn(M);
          
            % Notify me that it's now generating plots
            sprintf('Plotting...')
            
            % To generate figure files and forbid pop-ups
            % set(gcf, 'Visible', 'off'); 
            
            % Do a fourier transform and plot it
            % FFT
            % FT = fft(mag);
            % figure
            % FTplot
            
            % Plot the velocity as a function of time
            % figure
            % velocityt
            
            % Plot the phase space (theta vs alpha)
            % figure
            % phasespace
            
            % figure
            % threedplot

            % thetatime
            % recent
            % figure
            % alphaplot

            % figure
            % thetatime

            % figure
            % velocitymodtime
            
            % figure
            % sectionplot

            % figure
            % distanceplot
            
            % Play animation, this will close all other figures
            figure
            animation
            
            % print('-depsc', '-tiff', '-r300', '/home/matt/Dropbox/Billiards/Figures/velocityt')
            % print('-depsc', '-tiff', '-r300', 'C:/Users/Me/Dropbox/Billiards/Figures/velocityt')

                % Trying to generate figures for LaTeX
%  %              SAVE FIGURE LAPRINT
%  %              figfontsizes(11, 9);
%  %              laprint(gcf, 'fig', 'width', 9, ...
%  %              'asonscreen', 'off', 'keepfontprops', 'on', ...
%  %              'factor', 1, 'scalefonts', 'off', 'mathticklabels', 'on')
%  %              
%  %              Save Figures: Figures are labeled by 4 numbers, k_epsilon_w_vx
%  %              print(fig, '-djpeg', sprintf('/home/matt/Documents/BilliardProject/Dec2013/test1/Figure_% d_% d_% d_% d.eps', k, epsilon, thetarange, alpharange))
%  %             
%  %              tikssave = sprintf('/home/matt/Documents/BilliardProject/Dec2013/test1/Figure_% d_% d_% d_% d.tikz', k, epsilon, thetarange, alpharange)'
%  %              matlab2tikz('tikssave');
%  %              
%  %              Data Set: saved and labeled by 4 numbers, k_epsilon_w_vx
%  %              savefile = sprintf('/home/matt/Documents/BilliardProject/Dec2013/test1/File_% d_% d_% d_% d', k, epsilon, thetarange, alpharange);
%  %              save(savefile)      
        end
    end
    
    % Find average slope/linear fit for the given epsilon
    % avgslope(h) = mean(mean(slope(h, :, :)));
    % avglinslope(h) = mean(mean(linslope(h, :, :)));
    % lastvavg(h) = mean(mean(lastvarray(h, :, :)));
end

% Save a file with the average slopes for each epsilon
% save finalavgs lastvavg avglinslope avgslope


toc  % How long did the program take?