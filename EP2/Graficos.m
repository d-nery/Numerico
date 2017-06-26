%% Parte 1
X = dlmread('out1.txt');
t = X(:,1);
X = X(:,2);

xt = @(t) t + 1/(1 - t);

figure;
hold on;
plot(t, X, '.');
fplot(xt, [t(1) t(end)]);

title('Teste 1');
xlabel('t');
ylabel('x');

legend('Pontos calculados', 'Valor real');

%% Parte 2
X = dlmread('out2.txt');
t = X(:,1);
X = X(:,2:end);

xt = {@(t)  exp(-t)*sin(t) + exp(-3*t)*cos(3*t);
      @(t)  exp(-t)*cos(t) + exp(-3*t)*sin(3*t);
      @(t) -exp(-t)*sin(t) + exp(-3*t)*cos(3*t);
      @(t) -exp(-t)*cos(t) + exp(-3*t)*sin(3*t)};

figure;
hold on;
for j = 1:4
    plot(t, X(:,j), '.');
    fplot(xt{j}, [t(1) t(end)]);
end

title('Teste 2');
xlabel('t');
ylabel('x');

legend('X1 calculado', 'X1 real', 'X2 calculado', 'X2 real', 'X3 calculado', 'X3 real', 'X4 calculado', 'X4 real');

%% Parte 3
X = dlmread('out3.txt');
t = X(:,1);
X = X(:,2:end);

m = 7;

xt = {};
for i = 1:m
   xt{i} = @(t) exp(-(2*(1-cos(pi()/(m + 1))))*t)*sin(pi()*i/(m + 1))+exp(-(2*(1-cos(m*pi()/(m + 1))))*t)*sin(m*pi()*i/(m + 1));
end

figure;
hold on;
for j = 1:7
    plot(t, X(:,j), '.');
    fplot(xt{j}, [t(1) t(end)]);
end

title('Teste 3');
xlabel('t');
ylabel('x');

legend('X1 calculado', 'X1 real', ...
    'X2 calculado', 'X2 real', ...
    'X3 calculado', 'X3 real', ...
    'X4 calculado', 'X4 real', ...
    'X5 calculado', 'X5 real', ...
    'X6 calculado', 'X6 real', ...
    'X7 calculado', 'X7 real');

%% Parte 4 - Chua
f1 = figure;
s1 = subplot(4,3,[1 2 3 4 5 6 7 8 9]);
s2 = subplot(4,3,10);
s3 = subplot(4,3,11);
s4 = subplot(4,3,12);

for n = 1400:2499
    X = dlmread(['out4_' int2str(n) '.txt']);
    t   = X(:,1);
    Vc1 = X(:,2);
    Vc2 = X(:,3);
    Il  = X(:,4);
    
    plot3(s1, Vc1, Vc2, Il, '-b');
    view([45 45]);
    title(s1, ['R = ' int2str(n)]);
    xlabel(s1, 'VC1');
    ylabel(s1, 'VC2');
    zlabel(s1, 'IL');
    
    plot(s2, t, Vc1, '-b');
    ylabel(s2, 'VC1');
    xlabel(s2, 'Tempo (s)');

    
    plot(s3, t, Vc2, '-b');
    ylabel(s3, 'VC2');
    xlabel(s3, 'Tempo (s)');
    
    plot(s4, t, Il, '-b');
    ylabel(s4, 'IL');
    xlabel(s4, 'Tempo (s)');
    
    pause(0.1);
end

%% Tempos de execuçcao
T = dlmread('tempos.txt');
R = T(:,1);
T = T(:,2);
figure;
plot(R, T, '-');
title('Tempos de execução');
xlabel('R (ohm)');
ylabel('Tempo (s)');