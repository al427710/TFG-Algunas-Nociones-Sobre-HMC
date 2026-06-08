clear
clc

I_exacta = 1;

Nvals = [10 50 100 500 1000 5000 10000 50000 100000 500000];

M = 100; % numero de simulaciones

errores = zeros(size(Nvals));

for k = 1:length(Nvals)

    N = Nvals(k);

    error_medio = 0;

    for j = 1:M

        x = rand(N,1);

        Ihat = mean(x.*exp(x));

        error_medio = error_medio + abs(I_exacta - Ihat);

    end

    errores(k) = error_medio/M;

end

teorica = 1./sqrt(Nvals);

figure

loglog(Nvals,errores,'o-','LineWidth',2)
hold on

loglog(Nvals,teorica,'r--','LineWidth',2)

grid on

xlabel('Numero de muestras N')
ylabel('Error medio')

legend('Error Monte Carlo','1/sqrt(N)')

title('Convergencia Monte Carlo')
