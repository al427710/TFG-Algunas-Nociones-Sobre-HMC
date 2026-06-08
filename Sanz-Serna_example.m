clear
clc
close all

rng(1)

k = 100;
V = @(x) 0.5*k*(norm(x)-1)^2;
gradV = @(x) k*(norm(x)-1)*(x/norm(x));
N = 2000;
x0 = [1;0;0];

%% RANDOM WALK METROPOLIS

h = 0.5;
Xrw = zeros(3,N);
Xrw(:,1)=x0;
rejRW = 0;
for n=2:N
    x = Xrw(:,n-1);
    xprop = x + h*randn(3,1);
    alpha = exp(-(V(xprop)-V(x)));
    if rand < min(1,alpha)
        Xrw(:,n)=xprop;
    else
        Xrw(:,n)=x;
        rejRW = rejRW+1;
    end
end


%% HMC

T = 1;
dt = 0.1;
L = round(T/dt);
Xhmc = zeros(3,N);
Xhmc(:,1)=x0;

rejHMC = 0;

for n=2:N
    x = Xhmc(:,n-1);
    p = randn(3,1);
    xnew = x;
    pnew = p;

    % Leapfrog-Verlet
    pnew = pnew - 0.5*dt*gradV(xnew);
    for j=1:L
        xnew = xnew + dt*pnew;
        if j<L
            pnew = pnew - dt*gradV(xnew);
        end
    end
    pnew = pnew - 0.5*dt*gradV(xnew);

    H0 = 0.5*(p'*p) + V(x);
    H1 = 0.5*(pnew'*pnew) + V(xnew);
    alpha = exp(-(H1-H0));
    if rand < min(1,alpha)
        Xhmc(:,n)=xnew;
    else
        Xhmc(:,n)=x;
        rejHMC = rejHMC+1;

    end

end


%% AUTOCORRELACION



acfRW  = autocorr(Xrw(1,:),NumLags=20);
acfHMC = autocorr(Xhmc(1,:),NumLags=20);


%% FIGURA RW

figure
subplot(1,2,1)
plot(Xrw(1,:),Xrw(2,:),'.')
axis equal
xlabel('x_1')
ylabel('x_2')
title(sprintf('MH  (%d rechazos)',rejRW))
subplot(1,2,2)
stem(0:20,acfRW,'filled')
xlabel('Lag')
ylabel('Correlacion')
title('Autocorrelación MH')


%% FIGURA HMC


figure
subplot(1,2,1)
plot(Xhmc(1,:),Xhmc(2,:),'.')
axis equal
xlabel('x_1')
ylabel('x_2')
title(sprintf('HMC  (%d rechazos)',rejHMC))
subplot(1,2,2)
stem(0:20,acfHMC,'filled')
xlabel('Lag')
ylabel('Correlacion')
title('Autocorrelación HMC')

%% RESULTADOS


fprintf('\n');
fprintf('RANDOM WALK\n');
fprintf('Rechazos = %d\n',rejRW);
fprintf('Aceptación = %.2f %%\n',100*(1-rejRW/N));

fprintf('\n');

fprintf('HMC\n');
fprintf('Rechazos = %d\n',rejHMC);
fprintf('Aceptacion = %.2f %%\n',100*(1-rejHMC/N));