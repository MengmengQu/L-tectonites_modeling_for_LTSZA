% plot aspect ratio vs strain
% 2024.12.23

a1 = 1:0.1:100;
a2 = a1;
a3 = ones(size(a1));
r_a = a1./a3;

b1 = 1:0.1:100;
b2 = ones(size(b1));
b3 = 1./b1;
r_b = b1./b3;

Strain_a = sqrt((log(a1./a2)).^2+(log(a2./a3)).^2);
Strain_b = sqrt((log(b1./b2)).^2+(log(b2./b3)).^2);

plot(r_a,Strain_a,r_b,Strain_b)
xlabel('Aspect ratio');
ylabel('Strain');
xlim([1 100]);
axis square