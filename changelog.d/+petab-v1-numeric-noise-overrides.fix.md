PEtab v1: Fixed numeric `noiseParameters` overrides replacing the whole noise
formula instead of being substituted into it. This led to wrong sigmas and
likelihoods for noise formulas other than the bare noise placeholder
(#3560).
