clc;
clear;
close all;


%% GENERACIÓN DE DATOS
rng(1);
N = 1000;      % número de observaciones
d = 50;        % dimensión del problema
rho = 0.8;
SigmaX = zeros(d);
for i=1:d
    for j=1:d
        SigmaX(i,j) = rho^(abs(i-j));
    end
end
X = mvnrnd(zeros(d,1),SigmaX,N);
% Parámetros reales
beta_true = randn(d,1);
% Función sigmoide
sigmoid = @(z) 1 ./ (1 + exp(-z));
% Probabilidades reales
p_true = sigmoid(X*beta_true);
% Variables binarias
y = rand(N,1) < p_true;

sigma2 = 10;

%% ENERGÍA POTENCIAL Y GRADIENTE
V = @(beta) -sum( y .* log(sigmoid(X*beta) + 1e-12) + (1-y).*log(1 - sigmoid(X*beta) + 1e-12) ) + (1/(2*sigma2))*(beta'*beta);
gradV = @(beta) X'*(sigmoid(X*beta)-y) + beta/sigma2;

%% PARÁMETROS GENERALES
Nsamples = 5000;

%% METROPOLIS-HASTINGS

fprintf('Ejecutando Metropolis-Hastings...\n');
beta_MH = zeros(d,Nsamples);
step_MH = 0.03;
accept_MH = 0;

for i = 2:Nsamples
    current = beta_MH(:,i-1);
    proposal = current + step_MH*randn(d,1);
    log_alpha = -V(proposal) + V(current);
    if log(rand) < log_alpha
        beta_MH(:,i) = proposal;
        accept_MH = accept_MH + 1;
    else
        beta_MH(:,i) = current;
    end
end
acc_rate_MH = accept_MH / Nsamples;

%% HAMILTONIAN MONTE CARLO

fprintf('Ejecutando HMC...\n');
beta_HMC = zeros(d,Nsamples);
eps = 0.01;
L = 40;
accept_HMC = 0;

for i = 2:Nsamples

    q_current = beta_HMC(:,i-1);
    p_current = randn(d,1);
    q = q_current;
    p = p_current;
    %% LEAPFROG
    % Medio paso momento
    p = p - 0.5*eps*gradV(q);
    for j = 1:L
        % Paso completo posición
        q = q + eps*p;
        % Paso completo momento
        if j ~= L
            p = p - eps*gradV(q);
        end
    end
    % Último medio paso momento
    p = p - 0.5*eps*gradV(q);

    % Reversibilidad
    p = -p;
    
    %% HAMILTONIANOS
    current_H = V(q_current) + 0.5*(p_current'*p_current);
    proposed_H = V(q) + 0.5*(p'*p);
    log_alpha = current_H - proposed_H;

    %% ACEPTAR / RECHAZAR
    if log(rand) < log_alpha
        beta_HMC(:,i) = q;
        accept_HMC = accept_HMC + 1;
    else
        beta_HMC(:,i) = q_current;
    end
end
acc_rate_HMC = accept_HMC / Nsamples;


%% RESULTADOS NUMÉRICOS
fprintf('\n');
fprintf('=============================\n');
fprintf('RESULTADOS\n');
fprintf('=============================\n');
fprintf('MH Acceptance Rate  = %.3f\n',acc_rate_MH);
fprintf('HMC Acceptance Rate = %.3f\n',acc_rate_HMC);


%% TRACEPLOTS
figure;

subplot(2,1,1);
plot(beta_MH(1,:),'LineWidth',1);
title('Trayectoria Metropolis-Hastings');
xlabel('Iteración');
ylabel('\beta_1');

subplot(2,1,2);
plot(beta_HMC(1,:),'LineWidth',1);
title('Trayectoria HMC');
xlabel('Iteración');
ylabel('\beta_1');


%% AUTOCORRELACIÓN
figure;

subplot(2,1,1);
autocorr(beta_MH(1,:),NumLags=100);
title('Autocorrelación MH');

subplot(2,1,2);
autocorr(beta_HMC(1,:),NumLags=100);
title('Autocorrelación HMC');

%% SCATTER PLOTS
figure;

subplot(1,2,1);
scatter(beta_MH(1,:),beta_MH(2,:),5,'filled');
title('MH');

subplot(1,2,2);
scatter(beta_HMC(1,:),beta_HMC(2,:),5,'filled');
title('HMC');


%% EFFECTIVE SAMPLE SIZE (ESS)
maxLag = 100;

acf_MH = autocorr(beta_MH(1,:),NumLags=maxLag);
ess_MH = Nsamples / (1 + 2*sum(acf_MH(2:end)));
acf_HMC = autocorr(beta_HMC(1,:),NumLags=maxLag);
ess_HMC = Nsamples / (1 + 2*sum(acf_HMC(2:end)));

fprintf('\n');
fprintf('ESS MH  = %.2f\n',ess_MH);
fprintf('ESS HMC = %.2f\n',ess_HMC);


%% COMPARACIÓN FINAL
figure;
bar([ess_MH ess_HMC]);
set(gca,'XTickLabel',{'MH','HMC'});
ylabel('Effective Sample Size');
title('Comparación ESS');
grid on;