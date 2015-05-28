neighbour=zeros(1,M);
cmin=0;
cmax=2*pi;
range=0.1;

for n=1:M
    for j=1:M
        if (modt(j)-modt(n))^2 + (mag(j)-mag(n))^2 <= range^2
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

caxis([cmin, cmax]);
scatter(modt, mag, 3, neighbour);
axis([0 2*pi 0 vmax ]);
 

xlabel('Time mod(2pi)');
ylabel('Speed');