function f=collision(t,tn,xn,v, k,epsilon,omega)
% Tell if the particle is currently inside, or outside the billaird

x = xn + v*(t-tn);

if x(1) > 0 
    theta= atan(x(2)/x(1));
elseif x (1) < 0 
    theta = atan(x(2)/x(1)) + pi;
elseif ((x(1)==0) && (x(2)>0))
    theta= pi/2;
elseif ((x(1)==0) && (x(2)<0))
    theta=3*pi/2;
elseif ((x(1)==0) && (x(2)==0))
    theta=0;
end

%Positive Inside, Negitive Outside
f =  (1 + epsilon*sin(k*theta)*sin(omega*t)) - sqrt(x'*x);
