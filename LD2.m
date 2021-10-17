% Lukas Salavejus MNEfm-21
% Intelektualiuju sistemu 2LD

clear all; clc;

%20 skaiciu iejimo vektorius x
x = 0.1:1/22:1; 

%Duotoji funkcija
y = (1+0.6*sin(2*pi*x/0.7)+0.3*sin(2*pi*x))/2;
 plot(x,y);
 hold on;

%Struktura:
%1 pasleptasis sluoksnis
%5 neuronai

%Pirmojo (pasleptojo) sluoksnio pradiniai rysiu svoriai
w1 = [randn(1) randn(1) randn(1) randn(1) randn(1)];
b1 = [randn(1) randn(1) randn(1) randn(1) randn(1)];


%Antrojo sluoksnio (isejimo) pradiniai rysiu svoriai
w2 = [randn(1) randn(1) randn(1) randn(1) randn(1)];
b2 = randn(1);

%mokymo zingsnis
n = 0.1;
%Pirmo sluoksnio tinklo atsakas
v1 = zeros(1,length(w1));
y1 = zeros(1,length(w1));

for n_cycle = 1:1000000
    
    for i=1:length(x)
        for j=1:length(w1)
            %Skaiciuojamas tinklo atsakas
            %Pirmo sluoksnio neuronams
            v1(j) = x(i) * w1(j) + b1(j);
            %Aktyvcijos funkcijos pritaikymas. Naudojama sigmoidine funkcija
            %y = 1/(1+exp(-x))
            y1(j) = 1/(1+exp(-v1(j)));
        end

        %Tinklo atsakas antro sluoksnio neuronui
        v2 = sum(w2.*y1)+b2;
        %Aktyvacijos funkcijai naudojama tiesine fukncija y=1*x;
        y2 = 1 * v2;
        %klaida
        e=y(i)-y2;

        %Rysiu svoriu atnaujinimas
        %Klaidos gradientas isejimo sluoksnio neuronui
        delta2 = e;

            delta1 = zeros(1,length(w1));
            for j=1:length(w1)
                %Klaidos gradientai pasleptojo sluoksnio neuronams
                delta1(j) = y1(j) * (1- y1(j))*delta2 * w2(j);
                %Atnaujinami isejimo sluoksnio svoriai
                w2(j) = w2(j) + n * delta2 * y1(j);
                %Atnaujinami pasleptojo sluoksnio svoriai
                w1(j) = w1(j) + n*delta1(j)*x(i);
                b1(j) = b1(j) + n*delta1(j);
            end
        b2 = b2 + n * delta2;
    end
end

x_new = [0.1:1/22:1];
Y = zeros (1, length (x_new));

for i = 1:length(x_new);
    for j=1:length(w1)
            %Skaiciuojamas tinklo atsakas
            %Pirmo sluoksnio neuronams
            v1(j) = x_new(i) * w1(j) + b1(j);
            %Aktyvcijos funkcijos pritaikymas. Naudojama sigmoidine funkcija
            %y = 1/(1+exp(-x))
            y1(j) = 1/(1+exp(-v1(j)));
    end
    %Tinklo atsakas antro sluoksnio neuronui
    v2 = sum(w2.*y1)+b2;
    %Aktyvacijos funkcijai naudojama tiesine fukncija y=1*x;
    y2 = 1 * v2;
    Y(i) = y2;
end

plot(x_new,Y,'ro')









