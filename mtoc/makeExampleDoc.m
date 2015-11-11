opts.imageFormat = 'png';
for ie = 1:6
    save('ie.mat','ie','opts')
    wrap_path = fileparts(which('amiwrap.m'));
    opts.format = 'html';
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '.m']),opts)
    load('ie.mat')
    opts.format = 'html';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '_syms.m']),opts)
    load('ie.mat')
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '.m']),opts)
    load('ie.mat')
    opts.format = 'latex';
    wrap_path = fileparts(which('amiwrap.m'));
    publish(fullfile(wrap_path,'examples',['example_' num2str(ie)],['example_model_' num2str(ie) '_syms.m']),opts)
    load('ie.mat')
end
