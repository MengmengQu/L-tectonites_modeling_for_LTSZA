function [Strain, Flinnk, eigV_normalized] = Ltectonites_withOrientation(epsilon, gamma_t, gamma_dy, gamma_dx, t) 

% L tectonites simulation


L = [0 -gamma_t -gamma_dx; 0 epsilon gamma_dy; 0 0 -epsilon];
d =  0.5*(L + L');
d_n = norm(d)/sqrt(2);
L_n = L/d_n;

I = eye(3);
F = (I + L_n*(t/100))^100;   % From L to F

V = F*F';
[eigV, eigD] = eig(V);
[lambda, ind ] = sort(diag(eigD), 'descend');
a = sqrt(lambda);
eigV = eigV(:, ind);
eigV_normalized = eigV ./ sqrt(sum(eigV.^2, 1));

Strain = sqrt((log(a(1)/a(2)))^2+(log(a(2)/a(3)))^2);
Flinnk =  (log(a(1)/a(2))) / (log(a(2)/a(3)));

end



