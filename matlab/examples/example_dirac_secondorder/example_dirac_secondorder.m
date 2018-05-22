function example_dirac_secondorder()
%%
% COMPILATION

[exdir,~,~]=fileparts(which('example_dirac_secondorder.m'));
% compile the model
amiwrap('model_dirac_secondorder','model_dirac_secondorder_syms',exdir,1)

%%
% SIMULATION

% time vector
t = linspace(0,3,1001);
p = [1;0.5;2;3];
k = [];

options = amioption('sensi',0,...
    'maxsteps',1e4);
options.sens_ind = [3,1,2];

% load mex into memory
[msg] = which('simulate_model_secondorder_dirac'); % fix for inaccessability problems
options.sensi = 2;
sol = simulate_model_dirac_secondorder(t,log10(p),k,[],options);

%%
% FORWARD SENSITIVITY ANALYSIS

options.sensi = 2;

sol = simulate_model_dirac_secondorder(t,log10(p),k,[],options);

%%
% FINITE DIFFERENCES

optionsfd = options;
optionsfd.sensi = 1;
optionsfd.sens_ind = 1:4;

eps = 1e-4;
xi = log10(p);
solfd = simulate_model_dirac_secondorder(t,xi,k,[],optionsfd);
for ip = 1:length(options.sens_ind);
    xip = xi;
    xip(options.sens_ind(ip)) = xip(options.sens_ind(ip)) + eps;
    solp = simulate_model_dirac_secondorder(t,xip,k,[],optionsfd);
    s2x_fd(:,:,:,ip) = (solp.sx - solfd.sx)/eps;
    s2y_fd(:,:,:,ip) = (solp.sy - solfd.sy)/eps;
end

%%
% PLOTTING
if(usejava('jvm'))
    figure
    c_x = get(gca,'ColorOrder');
    for ip = 1:4
        for jp = 1:length(options.sens_ind)
            subplot(length(options.sens_ind),4,(jp-1)*4+ip)
            hold on
            for ix = 1:size(sol.x,2)
                plot(t,sol.s2x(:,ix,ip,jp),'.-','Color',c_x(ix,:))
                plot(t,s2x_fd(:,ix,ip,jp),'--','Color',c_x(ix,:))
            end
            ylim([-10,10])
            legend('x1','x1_{fd}','x2','x2_{fd}','Location','NorthEastOutside')
            legend boxoff
            title(['state sensitivity for p' num2str(ip) '-p' num2str(options.sens_ind(jp))])
            xlabel('time t')
            ylabel('x')
            box on
        end
    end
    set(gcf,'Position',[100 300 1200 500])
    figure
    for ip = 1:4
        for jp = 1:length(options.sens_ind)
            subplot(length(options.sens_ind),4,(jp-1)*4+ip)
            plot(t,abs(sol.s2x(:,:,ip,jp)-s2x_fd(:,:,ip,jp)),'r--')
            legend('error x1','error x2','Location','NorthEastOutside')
            legend boxoff
            title(['state sensitivity for p' num2str(ip) '-p' num2str(options.sens_ind(jp))])
            xlabel('time t')
            ylabel('error')
            ylim([1e-12,1e0])
            set(gca,'YScale','log')
            box on
        end
    end
    set(gcf,'Position',[100 300 1200 500])
    
    drawnow
end

end
