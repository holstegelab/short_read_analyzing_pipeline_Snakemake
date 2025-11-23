## Calculates population prior for variant based on genotype likelihoods (PL field)
## Simple version: includes genotype for which prior is calculated
## 'Non-simple' version: conditions on genotype for which prior is calculated.
##
## PROBLEM
## Genotype probability vectors cannot be used 'as-is'. 
## This can be best be explained from the case where no data is available: all genotypes are then equally likely, and therefore GATK gives a probability of 0.33/0.33/0.33. 
## One can observe two problems::
## - If this is a rare variant, then using a genotype probability of 0.33 here for an alternate homozygous allele is a large overestimation and would negatively affect the analysis significantly. 
## - Giving an equal probability to each genotype does not take into account Hardy-Weinberg.
## 
## 
## SOLUTION
## The functions in this file calculate genotype priors based on estimated pouplation allele frequency
## e.g. a low population frequency means also a lower prior probability for homozygous alternate genotypes. 
## It therefore estimates the allele frequency, takes this through hardy-weinberg, and then calculates a genotype prior. 
##  
## INPUT:
## matrix of (n_samples * n_genotypes) with genotype probabilities (based on PL field)
## - supported n_genotypes: 2 or 3
##
## OUTPUT:
## 1) allele frequency estimates
## 2) genotype priors
##    - simple: vector of n_genotypes with prior per genotype
##    - 'non-simple': vector of n_samples * n_genotypes with prior per genotype per sample
##   
## IMPLEMENTATION NOTES
##  1) Simple versus non-simple
##  A problem of this implementation is that it uses the genotype probabilities to estimate population frequencies, and then in turn uses the population frequency as prior for the genotype probabilities. 
##  This is circular. If n_samples is very high this is however not much of a problem, as each call only contributes a 1/n fraction to the genotype prior. Therefore, the prior for a specific call is for (n-1)/n part 
##  based on other calls (non-circular), and should therefore not be affected much by the specific call itself. 
##  However, for rare variants one could reason that there still might be an impact. The 'non-simple' version (which is not used in practice, as it has some other problems) therefore calculates for each call a separate prior, which does
##  not take into account the probability of the specific call itself at all (only the probability of the other calls). 
##
##  2) Chicken/egg problem
##  Note that the allele frequency is based on the genotype probabilities. Here we cannot use the direct genotype probabilities from GATK due to the abovementioned problem, as this would severely overstimate/underestimate 
##  allele-frequency in many cases. Allele frequency has to be estimated based on the posterior genotype probabilities. As these are not yet available (we need to calculate the population prior for that), this is somewhat 
##  of a problem. To solve this chicken/egg problem, we start with an allele frequency estimate of 1/n_allele, and calculate a posterior probability based on that. Then we refine our allele frequency estimate based 
##  on these posterior probabilities, and recalculate new psoterior genotype probabilities using this better estimate of the allele frequency. We iterate this until convergence. 


import numpy

def estimate_afs_simple(vectors, start=None):
    nsamples,ngeno = vectors.shape
    assert ngeno == 2 or ngeno == 3
    nallele = 2

    has_info = (vectors != (1.0 / ngeno)).all(axis=1)
    vectors = vectors[has_info,:]


    if start is None:
        afs_estimate = numpy.ones(nallele,dtype=float) / nallele
    else:
        afs_estimate = start
    genotype_priors = numpy.zeros(ngeno, dtype=float)
    max_diff = 1.0 


    if ngeno == 2:
        while(max_diff > 0.0000001):
            genotype_priors[0] = afs_estimate[0]
            genotype_priors[1] = afs_estimate[1]


            
            gvectors = genotype_priors *vectors
            gvectors = gvectors / numpy.sum(gvectors,axis=1)[:,numpy.newaxis]
            total_prob = numpy.sum(gvectors,axis=0)

            allele_counts = numpy.zeros(nallele,dtype=float)
            for j in range(nallele):
                idx = j
                allele_counts[j] += total_prob[idx]

            new_afs_estimate = allele_counts / float(len(vectors))

            max_diff = numpy.max(new_afs_estimate - afs_estimate)
            afs_estimate = new_afs_estimate

        genotype_priors[0] = afs_estimate[0]
        genotype_priors[1] = afs_estimate[1]
    else:
        while(max_diff > 0.0000001):
            for j in range(nallele):
                for k in range(j, nallele):
                    idx = int(0.5 * k * (k+1) + j)
                    genotype_priors[idx] = ((j!=k) + 1.) * afs_estimate[j] * afs_estimate[k]
            
            gvectors = genotype_priors *vectors
            gvectors = gvectors / numpy.sum(gvectors,axis=1)[:,numpy.newaxis]
            total_prob = numpy.sum(gvectors,axis=0)

            allele_counts = numpy.zeros(nallele,dtype=float)
            for j in range(nallele):
                for k in range(j, nallele):
                    idx = int(0.5 * k * (k+1) + j)
                    allele_counts[j] += total_prob[idx]
                    allele_counts[k] += total_prob[idx]

            new_afs_estimate = allele_counts / (len(vectors) * 2.0)

            max_diff = numpy.max(new_afs_estimate - afs_estimate)
            afs_estimate = new_afs_estimate

        for j in range(nallele):
            for k in range(j, nallele):
                idx = int(0.5 * k * (k+1) + j)
                genotype_priors[idx] = ((j!=k) + 1.) * afs_estimate[j] * afs_estimate[k]
        

    return (afs_estimate, genotype_priors[:,numpy.newaxis])


def estimate_afs(vectors):
    cache = {}
    nsamples,ngeno = vectors.shape
    nallele = 2 if ngeno == 3 else 1

    fil = numpy.ones(nsamples,dtype=bool)
    afs_estimates0 = numpy.zeros((nsamples, nallele),dtype=float)
    afs_estimates1 = numpy.zeros((nsamples, nallele),dtype=float)
    afs_estimates2 = numpy.zeros((nsamples, nallele),dtype=float)

    no_info = (vectors == (1.0 / ngeno)).all(axis=1)
    fil[no_info] = False
    
   
    v = vectors[fil,:]
    start_afs_estimate, xgpriors = estimate_afs_simple(v)

    start0 = start_afs_estimate
    start1 = start_afs_estimate
    start2 = start_afs_estimate

    for i in range(nsamples):
        plcol = tuple(vectors[i,:])
        if (i >= 1 and i < 10) or i ==  100:
            start0 = afs_estimates0[:i,:].mean(axis=0)
            start1 = afs_estimates1[:i,:].mean(axis=0)
            start2 = afs_estimates2[:i,:].mean(axis=0)
            
        if not plcol in cache:
            tfil = fil.copy()
            tfil[i] = True
            r = vectors.copy()
            r[i,:] = [1,0,0]
            afs_estimate0, gpriors = estimate_afs_simple(r[tfil,:], start=start0)
            r[i,:] = [0,1,0]
            afs_estimate1, gpriors = estimate_afs_simple(r[tfil,:], start=start1)
            r[i,:] = [0,0,1]
            afs_estimate2, gpriors = estimate_afs_simple(r[tfil,:], start=start2)

            cache[plcol] = (afs_estimate0, afs_estimate1, afs_estimate2)

        res = cache[plcol]
        afs_estimates0[i,:] = res[0]
        afs_estimates1[i,:] = res[1]
        afs_estimates2[i,:] = res[2]
        
    #estimate priors conditional per sample
    genotype_priors_per_sample = numpy.zeros((nsamples, int(0.5 * (nallele - 1) * nallele + nallele)), dtype=float)
   
    #3. determine genotype priors per sample
    for j in range(nallele):
        for k in range(j, nallele):
            idx = int(0.5 * k * (k+1) + j)
            if (j + k) == 0:
                caf = afs_estimates0
            elif (j + k) == 1:
                caf = afs_estimates1
            else:
                caf = afs_estimates2

            genotype_priors_per_sample[:,idx] = ((j!=k) + 1.) * caf[:,j] * caf[:,k]
    genotype_priors_per_sample = genotype_priors_per_sample / genotype_priors_per_sample.sum(axis=1)[:, numpy.newaxis]

    return (start_afs_estimate, genotype_priors_per_sample.T)

def estimate_afs_mix(vectors):
    if vectors.shape[1] == 2:
       afs = numpy.mean(vectors, axis=0)
       gafs = afs[:,numpy.newaxis]
       return (afs, gafs)
    nallele = 2 #FIXME

    afs_estimate = numpy.ones(nallele,dtype=float) / nallele
    genotype_priors = numpy.zeros(int(0.5 * (nallele - 1) * nallele + nallele), dtype=float)
    max_diff = 1.0 
    

    #estimate afs
    while(max_diff > 0.0000001):
        for j in range(nallele):
            for k in range(j, nallele):
                idx = int(0.5 * k * (k+1) + j)
                genotype_priors[idx] = ((j!=k) + 1.) * afs_estimate[j] * afs_estimate[k]
        
        gvectors = genotype_priors *vectors
        gvectors = gvectors / numpy.sum(gvectors,axis=1)[:,numpy.newaxis]
        total_prob = numpy.sum(gvectors,axis=0)

        allele_counts = numpy.zeros(nallele,dtype=float)
        for j in range(nallele):
            for k in range(j, nallele):
                idx = int(0.5 * k * (k+1) + j)
                allele_counts[j] += total_prob[idx]
                allele_counts[k] += total_prob[idx]

        new_afs_estimate = allele_counts / (len(vectors) * 2.0)

        max_diff = numpy.max(new_afs_estimate - afs_estimate)
        afs_estimate = new_afs_estimate
    
    
    #estimate posterior
    for j in range(nallele):
        for k in range(j, nallele):
            idx = int(0.5 * k * (k+1) + j)
            genotype_priors[idx] = ((j!=k) + 1.) * afs_estimate[j] * afs_estimate[k]
    
    gvectors = genotype_priors *vectors
    gvectors = gvectors / numpy.sum(gvectors,axis=1)[:,numpy.newaxis]

    
    #estimate priors conditional per sample
    genotype_priors_per_sample = numpy.zeros((len(vectors), int(0.5 * (nallele - 1) * nallele + nallele)), dtype=float)
   
    #1. determine allele counts per sample
    allele_counts_per_sample= numpy.zeros((len(vectors) , nallele),dtype=float)
    conditional_allele_freq1 = numpy.zeros((len(vectors), nallele),dtype=float)
    conditional_allele_freq2 = numpy.zeros((len(vectors), nallele),dtype=float)
    
    for j in range(nallele):
        for k in range(j, nallele):
            idx = int(0.5 * k * (k+1) + j)
            allele_counts_per_sample[:,j] += gvectors[:,idx]
            allele_counts_per_sample[:,k] += gvectors[:,idx]
   
    #2. determine conditional allele freq estimate
    for j in range(nallele):
        conditional_allele_freq1[:,j] = afs_estimate[j] + (1.0 - allele_counts_per_sample[:,j]) / (2.0 * len(vectors))
        conditional_allele_freq2[:,j] = afs_estimate[j] + (2.0 - allele_counts_per_sample[:,j]) / (2.0 * len(vectors))

    #3. determine genotype priors per sample
    for j in range(nallele):
        for k in range(j, nallele):
            idx = int(0.5 * k * (k+1) + j)
            if(j ==k):
                caf = conditional_allele_freq2
            else:
                caf = conditional_allele_freq1
            genotype_priors_per_sample[:,idx] = ((j!=k) + 1.) * caf[:,j] * caf[:,k]
        
    #gvectors = genotype_priors_per_sample *vectors
    #gvectors = gvectors / numpy.sum(gvectors,axis=1)[:,numpy.newaxis]


    return (afs_estimate, genotype_priors_per_sample.T)

def process_pl_ext(plrow, func=estimate_afs_simple):
    plrow = plrow - plrow.min(axis=0)
    plrow = 10.0**(plrow / -10.0) 
    plrow = plrow / plrow.sum(axis=0)
    afs, gafs = func(plrow.T)


    pdata = gafs * plrow
    posterior_cor = pdata / pdata.sum(axis=0)

    return (afs, gafs, posterior_cor)


def process_pl(plrow, func=estimate_afs_simple):
    plrow = plrow - plrow.min(axis=0)
    plrow = 10.0**(plrow / -10.0) 
    plrow = plrow / plrow.sum(axis=0)
    afs, gafs = func(plrow.T)


    pdata = gafs * plrow
    posterior_cor = pdata / pdata.sum(axis=0)

    return (afs, posterior_cor)
