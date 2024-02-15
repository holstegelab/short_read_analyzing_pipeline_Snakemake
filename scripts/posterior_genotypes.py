## Calculates population prior for variant based on genotype likelihoods (PL field)

import numpy
import cyvcf2
import sys
import numpy

#usage:  python get_posterior_dosage.py  input.bcf/vcf  chrom:from-start    output.bcf/vcf

def estimate_afs(vectors, ploidy=2, start=None):
    nsamples,ngeno = vectors.shape
    assert ploidy == 2 or ploidy == 1
    if ploidy == 2:
        nallele = numpy.cast[int]((-1 + numpy.sqrt(1 + 8 * ngeno)) / 2)
    else:
        nallele = ngeno

    has_info = (vectors != (1.0 / ngeno)).all(axis=1)
    vectors = vectors[has_info,:]


    #allele frequency vector
    if start is None:
        afs_estimate = numpy.ones(nallele,dtype=float) / nallele
    else:
        afs_estimate = start

    #genotyep frequency fector
    genotype_priors = numpy.zeros(ngeno, dtype=float)

    max_diff = 1.0 
    if ploidy == 1: #haploid
        while(max_diff > 0.0000001):
            for i in range(ngeno):
                genotype_priors[i] = afs_estimate[i]
            
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

        for i in range(ngeno):
            genotype_priors[i] = afs_estimate[i]
            
    else: #diploid
        while(max_diff > 0.0000001):
            #hardy weinberg
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



def process_pl(plrow, ploidy=2):
    plrow = plrow - plrow.min(axis=0)
    plrow = 10.0**(plrow / -10.0) 
    plrow = plrow / plrow.sum(axis=0)
    afs, gafs = estimate_afs(plrow.T, ploidy)


    pdata = gafs * plrow
    posterior_cor = pdata / pdata.sum(axis=0)

    return (afs, gafs, posterior_cor)



vcf_reader = cyvcf2.VCF(sys.argv[1])
vcf_reader.add_format_to_header({'ID':"GP", 'Number': 'G', 'Type':'Float', 'Description':'Probability (posterior)'})
vcf_reader.add_format_to_header({'ID':"PGQ", 'Number': '1', 'Type':'Float', 'Description':'GQ (posterior)'})
if len(sys.argv) > 2:
    out = cyvcf2.Writer(sys.argv[2], vcf_reader)
else:
    out = cyvcf2.Writer(sys.stdout, vcf_reader)

counter = 0
#walk through vcf
if len(sys.argv) > 3:
    threshold = int(sys.argv[3])
else:
    threshold = -1

maxlen = 1
for v in vcf_reader:
    counter +=  1
    if counter % 1000 == 0:
        sys.stderr.write(str(counter) + "\r")
        sys.stderr.flush()

    pls = v.format('PL')
    rad = v.format('AD')
    nallele = 1 + len(v.ALT)

    pl_missing = pls[:,0] <= -2147483647

    genotypes = v.genotypes
    ploidy = len(genotypes[0])-1
    missing_genotype = [-1] * ploidy + [False]


    pls[pl_missing,:] = 0 #set missing to 0
    rad[rad <= -2147483647] = 0 #set missing to 0
    
    depth = rad.sum(axis=1)
    depth[pl_missing] = 0

    if threshold == -1:
        afs, gafs, posterior_pls = process_pl(pls.T, ploidy=ploidy)
    else:
        fil = depth >= threshold
        afs, gafs, tposterior_pls = process_pl(pls[fil,:].T, ploidy=ploidy)
        
        pdata = gafs * pls.T
        posterior_pls = pdata / pdata.sum(axis=0)


    #convert to phred, with maximum value of 10000
    with numpy.errstate(divide='ignore'):
        phred_posterior_pls = -10.0 * numpy.log10(posterior_pls)
    phred_posterior_pls_min = phred_posterior_pls - phred_posterior_pls.min(axis=0)
    phred_posterior_pls_min[phred_posterior_pls_min > 10000.0] = 10000.0
    
    #get the best call
    calls = numpy.argmin(phred_posterior_pls_min,axis=0)
    
    #get the alleles for each genotype
    genotype_alleles = [None] * pls.shape[1]
    for j in range(nallele):
        for k in range(j, nallele):
            idx = int(0.5 * k * (k+1) + j)
            genotype_alleles[idx] = [j,k,False]

    #set the genotype per sample based on posterior calls
    for pos in range(len(calls)):
        if pl_missing[pos]:
            genotypes[pos] = missing_genotype
        else:
            genotypes[pos] = genotype_alleles[calls[pos]]
    v.genotypes = genotypes

    #reset missing to -2147483648
    phred_posterior_pls_min = numpy.cast[int](phred_posterior_pls_min)
    phred_posterior_pls_min[0,pl_missing] = -2147483648
    phred_posterior_pls_min[1:,pl_missing] = -2147483647


    #get the gq of the best call
    gq = numpy.sort(phred_posterior_pls_min,axis=0)[1,:]
    gq[gq > 99] = 99


    v.set_format('PGQ', gq)
    v.set_format('GP',posterior_pls.T.copy(order='C'))

    out.write_record(v)

sys.stderr.write("Processed " + str(counter) + " records\n")
sys.stderr.flush()

out.close();
vcf_reader.close();   
