%% Parte 1
X = dlmread('out1.txt');
t = X(:,1);
X = X(:,2);

xt = @(t) t + 1/(1 - t);

figure;
hold on;
plot(t, X, '.');
fplot(xt, [t(1) t(end)]);

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

%% Parte 4 - Chua
X = dlmread('out4.txt');
t   = X(:,1);
Vc1 = X(:,2);
Vc2 = X(:,3);
Il  = X(:,4);

figure;
hold on;
plot3(Vc1, Vc2, Il, '-b');
title('R = 1000');
xlabel('VC1');
ylabel('VC2');
zlabel('IL');

figure;
subplot(3,1,1);
hold on;
plot(t, Vc1, '-b');
title('R = 1000');
ylabel('VC1');

subplot(3,1,2);
plot(t, Vc2, '-b');
ylabel('VC2');

subplot(3,1,3);
plot(t, Il, '-b');
ylabel('IL');
xlabel('Tempo (s)');
