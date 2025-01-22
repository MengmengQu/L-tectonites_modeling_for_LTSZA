% test for L tectonites simulation
%      epsilon=1

% L = 0 -gamma_t -gamma_dx
%     0 epsilon  gamma_dy
%     0 0        -epsilon


clear

alpha    = 20;
t        = 0.1:0.1:6;
gamma_dy = 0:0.01:1;
gamma_dx = gamma_dy .* cotd(alpha);
gamma_t  = 0:0.1:10;
epsilon  = 1;


Strain  = zeros(length(t),length(gamma_dy),length(gamma_t));
Flinnk  = zeros(length(t),length(gamma_dy),length(gamma_t));
Orientation = zeros(3,3,length(t),length(gamma_dy),length(gamma_t));
phi = zeros(length(t),length(gamma_dy),length(gamma_t));
z_axis = [0; 0; 1];

for i = 1:length(t)
    for j = 1:length(gamma_dy)
        for k = 1:length(gamma_t)
            tt  = t(i);
            g_dy = gamma_dy(j);
            g_dx = gamma_dx(j);
            g_t = gamma_t(k);

           [Strain(i,j,k), Flinnk(i,j,k), Orientation(:,:,i,j,k)] = Ltectonites_withOrientation(epsilon, g_t, g_dy, g_dx, tt);
           PHI = rad2deg(acos(dot(Orientation(:,3,i,j,k),z_axis)));
           if PHI > 90
                PHI = 180 - PHI;
           end
           phi(i,j,k) = PHI;
        end
    end
end


% plot results    
gd_index   = find(gamma_dy == 0 | gamma_dy == 0.05 |  gamma_dy == 0.1|  gamma_dy == 0.2);
gamma_plot = [1,2,3,4,5];
plot_ori    = zeros(3,3,length(t),length(gd_index),length(gamma_plot));

figure
for jj = 1:length(gamma_plot)
    chosen_gamma_t = gamma_plot(jj);
    gt_index   = find(gamma_t == chosen_gamma_t);
    plotk = squeeze(Flinnk(:,gd_index,gt_index));
    plots = squeeze(Strain(:,gd_index,gt_index));

    subplot(2,3,jj)
    plot(plots(:,1),plotk(:,1),plots(:,2),plotk(:,2),plots(:,3),plotk(:,3),plots(:,4),plotk(:,4))
    xlabel('Strain');
    ylabel('K value');
    title('Plot of K value at \gamma_t = ', num2str(chosen_gamma_t));
    xlim([0 5])
end
legend('\gamma_d_y=0','\gamma_d_y=0.05','\gamma_d_y=0.1','\gamma_d_y=0.2')



figure
chosen_gamma_t_polar = 3;
gt_index_polar = find(gamma_t == chosen_gamma_t_polar);
strain_polar = squeeze(Strain(:,gd_index,gt_index_polar));  

for jj = 1:length(gd_index)

    plot_ori = squeeze(Orientation(:,:,:,gd_index(jj),gt_index_polar));
    [a1_evl, ~, a3_evl] = ConvertEigV2Angs(plot_ori);

    %     a1
      [~,a1in]  = find(a1_evl(2,:)<=(0.5*pi));
      [~,a1out] = find(a1_evl(2,:)>(0.5*pi));
      r1(a1in)  = sqrt(2) * sin(a1_evl(2,a1in)./2);
      r1(a1out) = sqrt(2) * cos(a1_evl(2,a1out)./2);


    %     a3
      [~,a3in]  = find(a3_evl(2,:)<(0.5*pi));
      [~,a3out] = find(a3_evl(2,:)>=(0.5*pi));
      r3(a3in)  = sqrt(2) * sin(a3_evl(2,a3in)./2);
      r3(a3out) = sqrt(2) * cos(a3_evl(2,a3out)./2);
      a3_evl(1,a3out) = a3_evl(1,a3out)-pi;


    subplot(2,2,jj)

      polarplot(a1_evl(1,a1in),r1(a1in),'.r',a1_evl(1,a1out),r1(a1out),'.r'...
          ,a1_evl(1,1),r1(1),'x',a1_evl(1,5),r1(5),'x',a1_evl(1,11),r1(11),'x',a1_evl(1,19),r1(19),'x',a1_evl(1,28),r1(28),'x'...
          ,a3_evl(1,a3in),r3(a3in),'.g',a3_evl(1,a3out),r3(a3out),'.g'...
          ,a3_evl(1,1),r3(1),'x',a3_evl(1,5),r3(5),'x',a3_evl(1,11),r3(11),'x',a3_evl(1,19),r3(19),'x',a3_evl(1,28),r3(28),'x')
    rlim([0 1])
    title('Plot of S_1 and S_3 evlutions at \gamma_d_y = ', num2str(gamma_dy(gd_index(jj))));

end

