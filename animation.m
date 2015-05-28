% Declare frames matrix
% No figures should be already open when movie generation 
close all

%  Initialization
fps = 100;                %  number of frames per second

figure
hold on;
axis equal; xlim([-2 2]); ylim = ([-2 2]);
pause on;
point = zeros(2,1);
temp = 0;

for n = (M-50):(M-1)
    time = tn(n+1) - tn(n);
    speed = 0.01;
    
    for j = 0:speed:time
        temp = temp + 1;
        theta = 0:.01:2*pi;
        b = polar(theta,1+epsilon.*sin(omega*(tn(n)+j)).*sin(k*theta),'-');
        title('Animation');
        % Parameterize the line
        point(:,1) = xn(:,n)+v(:,n)*(j);
        l = plot(point(1,1), point(2,1),'--r.','LineWidth',2,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','g',...
            'MarkerSize',15);
        pausetime = 0.0001;
        pause(pausetime);
        drawnow ;
        
        %  add frame to the movie
        FM(temp) = getframe(gcf);
        
        delete(b, l);
    end
           
end

% movie2avi(FM,'testanimate.avi','compression','none');


