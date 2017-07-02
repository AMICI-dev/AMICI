function example_neuron()
%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_neuron.m'));
% compile the model
amiwrap('model_neuron','model_neuron_syms',exdir)

%%
% SIMULATION

p = [0.02 0.3 65 0.9];

k = [-60,10];

% time vector
t = linspace(0,100,100);


options = amioption('sensi',0,...
    'maxsteps',1e5,...
    'nmaxevent',22);

% load mex into memory
[msg] = which('simulate_neuron_dirac'); % fix for inaccessability problems
sol = simulate_model_neuron(t,log10(p),k,[],options);

tic
sol = simulate_model_neuron(t,log10(p),k,[],options);

D = amidata(length(t),1,size(sol.z(sol.z<t(end)),2),size(sol.z(sol.z<t(end)),1),2);

rng(0);
D.Z = sol.z(sol.z<t(end));
% t = linspace(0,D.Z(end)-0.1,100);
D.Z = D.Z + 0.5*randn(size(D.Z));
D.Z(3) = NaN;
D.Sigma_Z = 0.5*ones(size(D.Z)); 
D.Z = D.Z + D.Sigma_Z.*randn(size(D.Z));

D.t = t;
D.condition = k;

sol = simulate_model_neuron([],log10(p),[],D,options);

llh = 0;
for ie = 1:size(D.Z,1)
    if(~isnan(D.Z(ie,1)))
        llh = llh - 1/2*log(2*pi*D.Sigma_Z(ie,1)^2) - 1/2*((sol.z(ie,1)-D.Z(ie,1))/D.Sigma_Z(ie,1))^2;
        if(sol.rz(ie,1)~=0)
            llh = llh - 1/2*log(2*pi*D.Sigma_Z(ie,1)^2) - 1/2*((sol.rz(ie,1))/D.Sigma_Z(ie,1))^2;
        end
    end
end

%%
% PLOTTING

if(usejava('jvm'))
    figure
    c_x = get(gca,'ColorOrder');
    subplot(2,1,1)
    for iz = 1:size(sol.x,2)
        plot(t,sol.x(:,iz),'.-','Color',c_x(iz,:))
        hold on
    end
    stem(sol.z,zeros(size(sol.z)))
    
    legend('x1','x2','events','Location','NorthEastOutside')
    legend boxoff
    xlabel('time t')
    ylabel('x')
    box on
    subplot(2,1,2)
    plot(-sol.llh,-llh,'bo')
    xlabel('amici nllh')
    ylabel('true nllh')
    set(gca,'XScale','log','YScale','log')
    hold on
    plot([1,1e4],[1,1e4],'k:')
    legend boxoff
    axis square
end


%%
% FORWARD SENSITIVITY ANALYSIS

options.sensi = 1;
options.sens_ind = 1:length(p);

sol = simulate_model_neuron(t,log10(p),k,D,options);

%%
% FINITE DIFFERENCES

eps = 1e-4;
xi = log10(p);
for ip = 1:4;
    xip = xi;
    xip(options.sens_ind(ip)) = xip(options.sens_ind(ip)) + eps;
    solp = simulate_model_neuron(t,xip,k,D,options);
    sz_fd(:,:,ip) = (solp.z - sol.z)/eps;
    sllh_fd(ip,1) = (solp.llh - sol.llh)/eps;
end

%%
% PLOTTING
if(usejava('jvm'))
    figure
    subplot(2,2,1)
    hold on
    plot(abs(sol.sz(:)),abs(sz_fd(:)),'bo')
    hold on
    plot([1e-5,1e5],[1e-5,1e5],'k:')
    xlim([1e-5,1e5])
    ylim([1e-5,1e5])
    title('event output sensitivity')
    xlabel('sz')
    ylabel('sz_{fd}')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    box on
    axis square
    
    subplot(2,2,2)
    plot(abs(sol.sz(:)),abs(sol.sz(:)-sz_fd(:)),'ro')
    hold on
    plot([1e-5,1e5],[1e-5,1e5],'k:')
    xlim([1e-5,1e5])
    ylim([1e-5,1e5])
    title('event output sensitivity')
    xlabel('sz')
    ylabel('error')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    box on
    axis square
    
    subplot(2,2,3)
    hold on
    plot(sol.sllh,sllh_fd,'ko')
    hold on
    plot([1e-5,1e5],[1e-5,1e5],'k:')
    xlim([1e-5,1e5])
    ylim([1e-5,1e5])
    xlabel('sllh')
    ylabel('sllh_{fd}')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    title('llh sensitivity')
    box on
    axis square
        
    subplot(2,2,4)
    plot(sol.sllh,abs(sol.sllh-sllh_fd),'ro')
    hold on
    plot([1e-5,1e5],[1e-5,1e5],'k:')
    xlim([1e-5,1e5])
    ylim([1e-5,1e5])
    xlabel('sllh')
    ylabel('error sllh')
    title('llh sensitivity')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    box on
    axis square
    set(gcf,'Position',[100 300 1200 500])
    
    drawnow
end

end