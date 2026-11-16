PEtab v2: Fixed numeric `noiseParameters` overrides replacing the whole noise
formula instead of being substituted into it. This led to wrong sigmas and
likelihoods for noise formulas other than the bare noise placeholder.
Timepoint-specific numeric overrides for such noise formulas now require
flattening (#3560).
