%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a scatter plot of PHASE SPACE
m=M;
neighbour=zeros(1,m);
cmin=0;
cmax=2*pi;
range=0.1;
% % 
for n=1:m
    for j=1:m
        if abs(thetan(j)-thetan(n))<=range && abs(alpha(j)-alpha(n))<=range
            neighbour(n)= neighbour(n)+1;
        end
    end
    
    if neighbour(n)>cmax
        cmax=neighbour(n);
    end
    if neighbour(n)<cmin
        cmin=neighbour(n);
    end
end
2
scatter(thetan, alpha,2, neighbour);
%plot(thetan, alpha,'.')

caxis([cmin, cmax]);
axis([0 2*pi 0 pi]); 
%title('Alpha and Theta');
xlabel('$\theta$ ');
ylabel('$\alpha$');