fit=p(2)*nn(:).^p(1);
linfit=ptwo(1).*nn(:) + ptwo(2);
% axis off

plot(nn, mag,'.', nn, linfit,'-')

xlabel('Number of Collisions');
ylabel('Speed');
% axis off