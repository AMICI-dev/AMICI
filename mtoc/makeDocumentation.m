setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin/'])
setenv('PATH', [getenv('PATH') ':/usr/local/bin'])
MatlabDocMaker.create('latex',true)

movefile(fullfile(wrap_path,'doc','latex','refman.pdf'),fullfile(wrap_path,'AMICI_guide.pdf'))
