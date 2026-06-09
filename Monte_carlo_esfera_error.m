clear
clc

rand('seed',12345);

% Volumen exacto de la esfera unidad
I_exacto = 4*pi/3;

% Valores de N a estudiar
Nvals = [100 500 1000 5000 10000 50000 100000 500000 1000000];

% Numero de repeticiones para cada N
M = 100;

errores = zeros(size(Nvals));

for k = 1:length(Nvals)

    N = Nvals(k);
    error_medio = 0;

    for j = 1:M

        % Puntos aleatorios en [-1,1]^3
        x = -1 + 2*rand(N,1);
        y = -1 + 2*rand(N,1);
        z = -1 + 2*rand(N,1);

        % Indicador de pertenencia a la esfera
        dentro = (x.^2 + y.^2 + z.^2 <= 1);

        % Estimacion Monte Carlo
        Ihat = 8*mean(dentro);

        % Acumulamos error absoluto
        error_medio = error_medio + abs(I_exacto - Ihat);

    end

    errores(k) = error_medio/M;

end

% Curva teorica O(N^{-1/2})
teorica = errores(1)*sqrt(Nvals(1))./sqrt(Nvals);

figure
loglog(Nvals,errores,'o-','LineWidth',2)
hold on
loglog(Nvals,teorica,'r--','LineWidth',2)

grid on
xlabel('Numero de muestras N')
ylabel('Error medio')
legend('Error Monte Carlo','O(N^{-1/2})','Location','southwest')
title('Convergencia volumen de la esfera')