clear
clc

rand('seed',12345);
I_exacto = 4*pi/3;

% PRIMERA APROXIMACIÓN UTILIZANDO 100 MUESTRAS
N1 = 100; 

% Generamos puntos aleatorios en el cubo [-1,1]^3
x = -1 + 2*rand(N1,1);
y = -1 + 2*rand(N1,1);
z = -1 + 2*rand(N1,1);

% Comprobamos cuantos puntos estan dentro de la esfera
dentro = (x.^2 + y.^2 + z.^2 <= 1);

% Numero de puntos dentro
puntos_dentro = sum(dentro);

% Volumen Montecarlo, volumen cubo = 2*2*2 = 8
Ihat1 = 8 * puntos_dentro / N1;
I1 = I_exacto-Ihat1

% SEGUNDA APROXIMACIÓN UTILIZANDO 50000 MUESTRAS
N2=50000;

x = -1 + 2*rand(N2,1);
y = -1 + 2*rand(N2,1);
z = -1 + 2*rand(N2,1);

dentro = (x.^2 + y.^2 + z.^2 <= 1);
puntos_dentro = sum(dentro);
Ihat2 = 8 * puntos_dentro / N2;

I2 = I_exacto-Ihat2