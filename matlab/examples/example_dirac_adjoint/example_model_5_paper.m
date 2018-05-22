clear
%% COMPILATION

[exdir,~,~]=fileparts(which('example_model_5.m'));
cd(exdir)
% compile the model
amiwrap('model_example_5','example_model_5_syms',exdir)

%% SIMULATION

% time vector
tout = linspace(0,4,9);
tfine = linspace(0,4,10001);
p = [1;0.5;2;3];
k = [];

D.Y = [  0.00714742903826096
      -0.00204966058299775
         0.382159034587845
          0.33298932672138
         0.226111476113441
         0.147028440865854
        0.0882468698791813
        0.0375887796628869
        0.0373422340295005];

D.Sigma_Y = 0.01*ones(size(D.Y));


options.sensi = 1;
options.sensi_meth = 'adjoint';
options.maxsteps = 2e4;
sol = simulate_model_example_5(tout,log10(p),k,D,options);
options.sensi = 0;
solfine = simulate_model_example_5(tfine,log10(p),k,[],options);

figure
errorbar(tout,D.Y,D.Sigma_Y)
hold on
plot(tfine,solfine.y)
legend('data','simulation')
xlabel('time t')
ylabel('observable')
title(['log-likelihood: ' num2str(sol.llh) ])

%% SYMBOLIC SOLUTION

syms x1(t) x2(t)
eqn1 = diff(x1) == -p(1)*x1;
eqn2 = x1(p(2)) == 1;
eqn3 = diff(x2) == p(3)*x1 - p(4)*x2;
eqn4 = x2(p(2)) == 0;

S = dsolve(eqn1,eqn2,eqn3,eqn4);
xx1 = matlabFunction(S.x1);
xx2 = matlabFunction(S.x2);

tfine = linspace(0,4,10001);
x1fine = xx1(tfine).*(tfine>=p(2));
x2fine = xx2(tfine).*(tfine>=p(2));

syms xB1(t) xB2(t)
eqn1B = diff(xB1) == p(1)*xB1 - p(3)*xB2;
eqn3B = diff(xB2) == p(4)*xB2;

syms sigma my 
x = sym('x',[2,1]);
J = -0.5*((x(2) - my)/sigma)^2;
dJdx = jacobian(J,x);
m0 = matlabFunction(dJdx);

mu1fine = zeros(size(x1fine));
mu2fine = zeros(size(x2fine));
for it = 1:length(tout);
    mm0 = m0(D.Y(it),D.Sigma_Y(it),x2fine(tfine == tout(it)));
    eqn2B = xB1(tout(it)) == mm0(1);
    eqn4B = xB2(tout(it)) == mm0(2);
    S = dsolve(eqn1B,eqn2B,eqn3B,eqn4B);
    xxB1 = matlabFunction(S.xB1);
    xxB2 = matlabFunction(S.xB2);
    mu1fine = mu1fine + xxB1(tfine).*(tfine<=tout(it));
    mu2fine = mu2fine + xxB2(tfine).*(tfine<=tout(it));
end

grad(1,1) = trapz(tfine,-x1fine.*mu1fine)*p(1)*log(10);
t_idx = find(tfine == p(2));
grad(2,1) = -(p(3)*mu2fine(t_idx)-p(1)*mu1fine(t_idx))*p(2)*log(10);
grad(3,1) = trapz(tfine,x1fine.*mu2fine)*p(3)*log(10);
grad(4,1) = trapz(tfine,-x2fine.*mu2fine)*p(4)*log(10);

%% FD

eps = 1e-7;
xi = log10(p);
grad_fd_f = NaN(4,1);
grad_fd_b = NaN(4,1);
for ip = 1:4;
    options.sensi = 0;
    xip = xi;
    xip(ip) = xip(ip) + eps;
    solpf = simulate_model_example_5(tout,xip,k,D,options);
    grad_fd_f(ip,1) = (solpf.llh-sol.llh)/eps;
    xip = xi;
    xip(ip) = xip(ip) - eps;
    solpb = simulate_model_example_5(tout,xip,k,D,options);
    grad_fd_b(ip,1) = -(solpb.llh-sol.llh)/eps;
end

figure
scatter(abs(grad_fd_f),abs(sol.sllh))
hold on
scatter(abs(grad_fd_b),abs(sol.sllh))
set(gca,'XScale','log')
set(gca,'YScale','log')
hold on
plot([1e1,1e4],[1e1,1e4],'k:')
xlim([1e1,1e4])
ylim([1e1,1e4])
legend('forward FD','backward FD')

%% figures paper

fs = 20;
lw = 2;
lwf = 1.5;
cm = lines;

% system

figure
plot(tout,D.Y,'x','LineWidth',2*lw,'Color',cm(5,:),'MarkerSize',15)
hold on
plot(tfine,solfine.y,'k','LineWidth',lw)
xlabel('time t')
ylabel('observable')
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
xlim([-0.1,4.1])
ylim([-0.05,0.4])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['y'])

% finite differences
figure
plot(tfine,solfine.y,'k:','LineWidth',lw)
hold on
xip = xi;
xip(1) = xip(1) + 1e-1;
solp = simulate_model_example_5(tfine,xip,k,[],options);
plot(tfine,solp.y,'LineWidth',lw,'Color',cm(1,:))
xlabel('time')
ylabel('observable')
xlim([-0.1,4.1])
ylim([-0.05,0.4])
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sy_fd'])


figure
bar(1:4,grad_fd_b,'LineWidth',lw,'FaceColor',[0.5,0.5,0.5])
set(gca,'XTick',[])
set(gca,'YTick',[])
ylabel('sensitivity J')
xlabel('parameters')
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
ylim([-300,300])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sJ_fd'])

% forward sensitivities
options.sensi = 1;
options.sensi_meth = 'forward';
options.maxsteps = 1e6;
% load mex into memory
sol = simulate_model_example_5(tfine,log10(p),k,[],options);

figure
hold on
plot(tfine,sol.sy(:,:,1),'LineWidth',lw,'Color',cm(2,:))
box on
xlabel('time')
ylabel('obs. sensitivity')
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
set(gca,'XTick',[])
set(gca,'YTick',[])
ylim([-0.6,0.4])
xlim([-0.1,4.1])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sy_fsa'])

sol = simulate_model_example_5(tout,log10(p),k,[],options);
sJ = -(squeeze(sol.sy)')*((sol.y-D.Y)./(D.Sigma_Y.^2));
figure
bar(1:4,sJ,'LineWidth',lw,'FaceColor',[0.5,0.5,0.5])
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('parameters')
ylabel('sensitivity J')
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
ylim([-300,300])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sJ_fsa'])

% adjoint sensitivities

figure
plot(tfine,mu1fine,'LineWidth',lw,'Color',cm(3,:))
hold on
plot(tfine,mu2fine,'LineWidth',lw,'Color',cm(3,:))
xlabel('time')
ylabel('adjoint state')
xlim([-0.1,4.1])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sy_asa'])

figure
bar(0,sum(-x1fine.*mu1fine*(4/(length(tfine)-1))*p(1)*log(10)),'LineWidth',lw,'FaceColor',[0.5,0.5,0.5])
hold on
plot(tfine,fliplr(cumsum(fliplr(-x1fine.*mu1fine*(4/(length(tfine)-1))*p(1)*log(10)))),'LineWidth',lw,'Color',cm(3,:))
scatter(0,sum(-x1fine.*mu1fine*(4/(length(tfine)-1))*p(1)*log(10)),40,cm(3,:),'o','LineWidth',lw)
xlabel('time')
ylabel('sensitivity J')
set(gca,'FontSize',fs)
set(gca,'LineWidth',lwf)
set(gca,'XTick',[])
set(gca,'YTick',[])
ylim([-300,300])
xlim([-0.1,4.1])
set(gcf,'PaperPositionMode','auto','Position',[100 300 300 200])
print('-depsc','-r300',['sJ_asa'])






