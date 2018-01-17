% ModelProblem
%
% Overview:
%   We calculate the resonance freqency for a small bubble with Dirichlet
%   boundary conditions on a bigger bubble, and compare against the
%   Minnaert resonance for the small bubble in free-space.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Matias Ruiz, Sanghyeon Yu, Erik Orvehed Hiltunen.

clear all
close all

distTol = 5e-6;
fTol = 1e-10;
iterMax = 100;
nPoints = 2^6;

kappa_w = 2*10^9;
rho_w = 1*10^3;

kappa_b = 140*10^3;
rho_b = 1.2;

v_w = sqrt(kappa_w/rho_w);
v_b = sqrt(kappa_b/rho_b);

delta = rho_b/rho_w;

radius_d = 0.2;
radius_b = 1;

initialGuess = 0.0001;
z0 = initialGuess;
z1 = initialGuess-initialGuess/100;
z2 = initialGuess-initialGuess/200;    

%%% compute the characteristic value of 'A' by applying Muller's method
M = 100;
r = linspace(0.03,1,M);
char_nums = zeros(M,1);
char_nums_A = zeros(M,1);
char_nums_freeSpace = zeros(M,1);
% char_nums_t = zeros(M,1);
% char_num_t = tools.MullersMethod(@(z) besselj(0,z), z0, z1, z2, iterMax, distTol, fTol); 
for n = 1:M    
    func2 = @(omega) g(omega, v_w, v_b, delta, radius_b, r(n), 0);
    char_nums(n) = tools.MullersMethod(func2, z0, z1, z2, iterMax, distTol, fTol);
    func = @(omega) f(omega, v_w, v_b, delta, radius_b, r(n), nPoints);
    char_nums_A(n) = tools.MullersMethod(func, char_nums(n)-0.0001, char_nums(n)+0.0001, char_nums(n), iterMax, distTol, fTol);
    func3 = @(omega) h(omega, v_w, v_b, delta, r(n), nPoints);
    char_nums_freeSpace(n) = tools.MullersMethod(func3, z0, z1, z2, iterMax, distTol, fTol);
%     char_nums_t(n) = char_num_t/(r(n)*v_b); 
end

% Plots the solution and saves good images
lfsz = 10;
fsz = 14;
alw = 0.75;
lw = 1.5;

plot(r(1:end-1),abs(real(char_nums(1:end-1))),r(1:end-1),abs(real(char_nums_freeSpace(1:end-1))))
xlabel('Radius R_d');
ylabel('Resonance frequency \omega_c');
hold on
% plot(r,real(char_nums_t))
set(gca, 'FontSize', fsz, 'LineWidth', alw);
l = legend('Model problem','Small bubble free-space', 'Small bubble Helmholtz eigenvalue','Location','East');
l.FontSize = lfsz;
hold off
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
print('simpleSol','-dpng');

figure
plot(r(8:end-6),abs(abs(real(char_nums(8:end-6)))-abs(real(char_nums_A(8:end-6))))./abs(real(char_nums(8:end-6))))
xlabel('Radius R_d');
ylabel('Relative error');set(gca, 'FontSize', fsz, 'LineWidth', alw);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
print('simpleErr','-dpng');