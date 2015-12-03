close all
clear
amiwrap('model_example_8','example_model_8_syms',pwd)

a = 2;
b = 1;
c = 1;
d = 1;

p = [a b c d];

v0 = 1;
I0 = 1;
tI = 1;

k = [v0,I0];

t = linspace(0,2.6,10002);

options.sensi = 1;

sol = simulate_model_example_8(t,p,k,[],options);

[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(p,@(theta) simulate_model_example_8(t,theta,k,[],options),1e-8,'x','sx');

figure
plot(t,sol.y)


figure
for ip = 1:length(p)
    subplot(length(p),1,ip)
    plot(t,g(:,:,ip),'.-')
    hold on
    plot(t,g_fd_f(:,:,ip),'.-')
    plot(t,g_fd_b(:,:,ip),'.-')
    plot(t,g_fd_c(:,:,ip),'.-')
    if(ip == 1)
        legend('g','fd_f','fd_b','fd_c','Location','EastOutside')
    end
end