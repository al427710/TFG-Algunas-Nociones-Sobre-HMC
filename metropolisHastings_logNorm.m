clc;
clear;
close all;

mu = 0;
sigma_ln = 1;
dist_obj = @(x) (1./x) .* exp(-(log(x)-mu).^2 ./ (2*sigma_ln^2)) .* (x > 0);

N = 50000;
sigma_q = 0.5;

x = zeros(N,1);
x(1) = 1;
aceptados = 0;

% Algoritmo Metropolis-Hastings
for i = 2:N
    x_current = x(i-1);
    x_prop = lognrnd(log(x_current), sigma_q);
    q_f = lognpdf(x_prop, log(x_current), sigma_q);
    q_b = lognpdf(x_current, log(x_prop), sigma_q);
    alpha = min(1, (dist_obj(x_prop) * q_b) / (dist_obj(x_current) * q_f));
    if rand < alpha
        x(i) = x_prop;
        aceptados = aceptados + 1;
    else
        x(i) = x_current;
    end
end

% Tasa de aceptación
tasa_acept = aceptados/(N-1);
fprintf('Tasa de aceptación: %.4f\n', tasa_acept);

% Burn-in
burn_in = 5000;
samples = x(burn_in:end);

% Histograma
figure;
histogram(samples, 100, 'Normalization', 'pdf');
hold on;

% Distribución real
xx = linspace(0.001,15,1000);
yy = lognpdf(xx,mu,sigma_ln);
plot(xx,yy,'r','LineWidth',2);
xlabel('x');
ylabel('Densidad');
title('Metropolis-Hastings - Distribución Log-Normal');
legend('Muestras MCMC','Distribución real');
grid on;

% Traza
figure;
plot(x);
xlabel('Iteración');
ylabel('Estado');
title('Trayectoria Cadena de Markov');
grid on;

% Estadísticas
fprintf('\n');
fprintf('Media estimada: %.4f\n', mean(samples));
fprintf('Varianza estimada: %.4f\n', var(samples));

% Valores teóricos
media_teorica = exp(mu + sigma_ln^2/2);
var_teorica = (exp(sigma_ln^2)-1) * exp(2*mu + sigma_ln^2);
fprintf('\n');
fprintf('Media teórica: %.4f\n', media_teorica);
fprintf('Varianza teórica: %.4f\n', var_teorica);
