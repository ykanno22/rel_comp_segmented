close all
clear
%
domain_x = 2.5*10^(0);
const.std = 0.25*10^(0);
const.f_per_x_elas = 15.0;
const.f_per_x_minus = 0.15 * const.f_per_x_elas;
const.f_per_x_plus  = 0.25 * const.f_per_x_elas;
num.division = 199;

list_of_x = (-domain_x:((2*domain_x)/num.division):domain_x)';
num.sample    = length(list_of_x);

rng(111,'twister');
w = domain_x;
list_of_noisy_x = -w + (2 * w * rand(num.division+1 ,1));
list_of_noisy_x = sort(list_of_noisy_x);
clear w
br_point(1) =     round(num.sample / 3);
br_point(2) = 2 * round(num.sample / 3);
rng(222,'twister');
r = const.std * randn(num.sample,1);
%
list_of_noisy_f = zeros(num.sample,1);
for ii=br_point(1):br_point(2)
    val_x = list_of_noisy_x(ii);
    list_of_noisy_f(ii) =...
        const.f_per_x_elas * val_x;
end
br_y(1) = min(list_of_noisy_f);
br_y(2) = max(list_of_noisy_f);
for ii=1:br_point(1)
    val_x = list_of_noisy_x(ii) - list_of_noisy_x(br_point(1));
    list_of_noisy_f(ii) =...
        br_y(1) + (const.f_per_x_minus * val_x);
end
for ii=br_point(2):num.sample
    val_x = list_of_noisy_x(ii) - list_of_noisy_x(br_point(2));
    list_of_noisy_f(ii) =...
        br_y(2) + (const.f_per_x_plus * val_x);
end
list_of_noisy_f = list_of_noisy_f + r;


list_of_noisy_f = 0.1 * list_of_noisy_f;

%%%% Scaling-->
list_of_noisy_x = 2.0 * list_of_noisy_x;
list_of_noisy_f = 2.0 * list_of_noisy_f;
%%%% <--Scaling

plot(list_of_noisy_x, list_of_noisy_f, 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3);
hold on;
grid on;
axis equal;
xlabel('Strain ($10^{-3}$ m/m)', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);

save('tri-modulus_1d_data_set.mat',...
    'list_of_noisy_x', 'list_of_noisy_f');

