%  Exercice 1b - Différences finies
%
%  Modélisation de la température en régime permanent
%
%  B. Giroux
%  INRS-ETE
%
%  2015-01-11
%

% température en surface
T0 = 280;   % K
% flux dans la N^e couche
Q = -55e-3;  % W/m^2

% nombre de couches
N = 20;
% épaisseur des couches
dz = 100;    % m

% Production de chaleur
A = zeros(N,1);
A(end-3:end) = 2.8e-6;  % W/m^3

% conductivité thermique  - W/m/K
l = zeros(N,1);
l(1:5) = 1.8;
l(6:10) = 3.7;
l(11:16) = 2.4;
l(end-3:end) = 3.5;


figure(1)
subplot(121)
plot(kron(l,[1;1]),[0;kron(((1:N-1)*dz)',[1;1]);N*dz])
set(gca,'YDir','reverse','XLim',[1.5 4])
xlabel('\lambda (W/m/K)','FontSize',14)
ylabel('Profondeur (m)','FontSize',14)

%%
%
% construction de la matrice A
%
M = zeros();
M(1,1) = 1;
for n=2:N
    M(n,n-1) =
    M(n,n) =
    M(n,n+1) =
end
M(N+1,N) = l(N)/dz;
M(N+1,N+1) = -l(N)/dz;

%
% vecteur source
%
b = zeros(N+1,1);
b(1) =
for n=2:N
    b(n) = -(A(n-1)+A(n))/2;
end
b(end) =

%
% Solution du système
%
x = M\b;

% vecteur des profondeurs
z=((0:N)*dz)';

figure(2)
subplot(122)
plot(x,z)
set(gca,'YDir','reverse')
xlabel('T (K)','FontSize',14)

%
% Affichage de la matrice
%
figure(3)
spy(M)

m = unique(M);
MM = M;
nn = 1:length(m);
label = cell(length(m),1);
for n=1:length(m)
    MM(M==m(n)) = n;
    label{n} = num2str(m(n));
end

figure(4)
imagesc(MM)
h=colorbar;
set(h, 'YTick', nn)
set(h,'YTickLabel', label)
axis equal
axis tight

figure(5)
plot(z,b-M*x,'o-')
title('Erreur = b - Ax','FontSize',18)
xlabel('Profondeur (m)','FontSize',14)
ylabel('Erreur','FontSize',14)
