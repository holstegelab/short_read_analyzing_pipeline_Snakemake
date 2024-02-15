import argparse

import cyvcf2
import numpy
import sys

def parse_gtw(gtw_sample):
    """Parse GTW field written by merge_phasing"""
    #format of gtw: pos + ref_len + allele1 + delimiter + allele2      
    pos = int(gtw_sample[:3])
    ref_len = 1
    if gtw_sample[3] == '+':
        #parse multiple digit prefix of gtw_sample[4:]
        #first determine how many digits there are
        ndigits = 1
        while gtw_sample[4+ndigits].isdigit():
            ndigits += 1
        ref_len = int(gtw_sample[4:4+ndigits])
        alleles = gtw_sample[4+ndigits:]
    else:
        alleles = gtw_sample[3:]
    if '!' in alleles:    
        alleles = alleles.split('!')
        delimiter = '!'
    else:
        alleles = alleles.split('|')
        delimiter = '|'

    return (pos, ref_len, alleles, delimiter)


def complete_pos(pos3, pos):
    """Positions are stored in the VCF using only the last 3 digits. This function completes the position
    by adding the prefix of the position to the last 3 digits. It returns the completed position and a boolean
    indicating whether the position is matching."""

    if pos3 == int(str(pos)[-3:]):
        return (pos, True)

    spos = str(pos)
    prefix_pos = int(spos[:-3]) * 1000
    
    #two possible positions
    x1 = prefix_pos + pos3 
    x2 = (prefix_pos - 1000 + pos3) 

    d1 = abs(x1 - pos)
    d2 = abs(x2 - pos)

    if d1 < d2:
        return (x1, False)
    else:
        return (x2, False)


def check_common_suffix(allele1_vcf, allele2_vcf, allele1_phased, allele2_phased, ref_len_phased, ref_len_vcf, stats):
    orientation= 0
    phase = False


    #print(f"DEBUG: {allele1_vcf}, {allele2_vcf}, {allele1_phased}, {allele2_phased}, {ref_len_phased}, {ref_len_vcf}")

    if ref_len_phased > ref_len_vcf: #if the phased reference is longer than the vcf reference, then we need to drop the common suffix from the phased alleles
        common_suffix = ref_len_phased - ref_len_vcf
        allele1_suffix = allele1_phased[-common_suffix:]
        allele2_suffix = allele2_phased[-common_suffix:]
        
        if allele1_suffix != allele2_suffix: #if the suffixes are not the same, then we cannot transfer phase
            print(f"DEBUG phased suffix mismatch: {allele1_vcf}, {allele2_vcf}, {allele1_phased}, {allele2_phased}, {ref_len_phased}, {ref_len_vcf}")
            return (orientation, phase)

        allele1_phased = allele1_phased[:-common_suffix]
        allele2_phased = allele2_phased[:-common_suffix]
        suffix_type = 'phased'


    elif ref_len_phased < ref_len_vcf: #if the phased reference is shorter than the vcf reference, then we need to drop the common suffix from the vcf alleles
        common_suffix = ref_len_vcf - ref_len_phased
        allele1_suffix = allele1_vcf[-common_suffix:]
        allele2_suffix = allele2_vcf[-common_suffix:]

        if allele1_suffix != allele2_suffix: #if the suffixes are not the same, then we cannot transfer phase
            print(f"DEBUG vcf suffix mismatch: {allele1_vcf}, {allele2_vcf}, {allele1_phased}, {allele2_phased}, {ref_len_phased}, {ref_len_vcf}")
            return (orientation, phase)

        allele1_vcf = allele1_vcf[:-common_suffix]
        allele2_vcf = allele2_vcf[:-common_suffix]

        suffix_type = 'vcf'

    #now check if the alleles match
    if allele1_phased == allele1_vcf and allele2_phased == allele2_vcf:
        phase = True
        orientation = 1
        stats[f'phase_transferred_after_dropping_{suffix_type}_suffix']  = stats.get(f'phase_transferred_after_dropping_{suffix_type}_suffix', 0) + 1
    elif allele1_phased == allele2_vcf and allele2_phased == allele1_vcf:
        phase = True
        orientation = 2
        stats[f'phase_transferred_after_dropping_{suffix_type}_suffix']  = stats.get(f'phase_transferred_after_dropping_{suffix_type}_suffix', 0) + 1



    #print("DEBUG RESULT: ", phase, orientation)
    return (orientation, phase)



def check_vcf_deletion(vcf_ref, allele1_vcf, allele2_vcf, allele1_phased, allele2_phased, preflen, vcf_pos, phased_pos, stats):
    orientation = 0
    phase = False


    #check if the phased variant covers the '*' allele in the VCF
    if preflen < (vcf_pos - phased_pos):
        return (orientation, phase)

    #reconstruct the phased alleles
    pseq1 = '*' * preflen
    pseq2 = '*' * preflen

    pseq1 = allele1_phased + pseq1[len(allele1_phased):]
    pseq2 = allele2_phased + pseq2[len(allele2_phased):]

    #now check if the alleles match
    diffpos = vcf_pos - phased_pos
    pseq1 = pseq1[diffpos:]
    pseq2 = pseq2[diffpos:]


    aseq1 = '*' * len(vcf_ref)
    aseq2 = '*' * len(vcf_ref)

    if allele1_vcf != '*':
        aseq1 = allele1_vcf + aseq1[len(allele1_vcf):]
 
    if allele2_vcf != '*':
        aseq2  = allele2_vcf + aseq2[:len(allele2_vcf):]



    min_len1 = min(len(pseq1), len(aseq1))
    min_len2 = min(len(pseq2), len(aseq2))

    if pseq1[:min_len1] == aseq1[:min_len1] and pseq2[:min_len2] == aseq2[:min_len2]:
        phase = True
        orientation = 1
        stats['phase_transferred_to_vcf_deletion']  = stats.get('phase_transferred_to_vcf_deletion', 0) + 1
    elif pseq1[:min_len2] == aseq2[:min_len2] and pseq2[:min_len1] == aseq1[:min_len1]:
        phase = True
        orientation = 2
        stats['phase_transferred_to_vcf_deletion']  = stats.get('phase_transferred_to_vcf_deletion', 0) + 1
    #if vcf_pos == 1014737:
    #    print("DEBUG INPUT: ", allele1_vcf, allele2_vcf, allele1_phased, allele2_phased, preflen, vcf_pos, phased_pos)
    #    print(f"DEBUG: {pseq1}, {pseq2}, {aseq1}, {aseq2}")
    #    util.debug_here()

    return (orientation, phase)


def check_double_deletion(vcf_ref, allele1_vcf, allele2_vcf, allele1_phased, allele2_phased, preflen, vcf_pos, phased_pos, stats):
    start_vcf = vcf_pos
    end_vcf = vcf_pos + len(vcf_ref)
    start_phased = phased_pos
    end_phased = phased_pos + preflen
    #check overlap
    if start_vcf > end_phased or start_phased > end_vcf:
        overlap = False
    else:
        overlap = True

    phase = False
    orientation = 0

    if overlap:
        #align * alleles
        if (allele1_vcf == '*' and allele1_phased == '*') or (allele2_vcf == '*' and allele2_phased == '*'):
            orientation = 1
            phase = True
        elif (allele2_vcf == '*' and allele1_phased == '*') or (allele1_vcf == '*' and allele2_phased == '*'):
            orientation = 2
            phase = True
        if phase:
            stats['phase_transferred_to_double_deletion']  = stats.get('phase_transferred_to_double_deletion', 0) + 1
    return (orientation, phase)

def process_record(vrecord, samples, stats, write_stack, last_failed):

    #first set all genotypes to unphased if they are phased
    idx_phased = numpy.where(vrecord.gt_phases)[0]

    if len(idx_phased) > 0:
        genotypes = vrecord.genotypes
        stats['nvariant_with_phased_gt'] = stats.get('nvariant_with_phased_gt', 0) + 1
        for i in idx_phased:
            a,b,p = vrecord.genotypes[i]
            genotypes[i] = [a,b,False]
        vrecord.genotypes = genotypes

    #short circuit if no GTW field
    if 'GTW' not in vrecord.FORMAT:
        return vrecord

    #walk across samples that have phasing information in gtw
    #get gtw / psw data
    stats['nvariant_with_gtw']  = stats.get('nvariant_with_gtw', 0) + 1
    

    gtw = vrecord.format('GTW')
    idx = numpy.where(gtw != '.')[0]

    if not len(idx) > 0: #no samples with phasing information
        return vrecord


    genotypes = vrecord.genotypes

    psw = vrecord.format('PSW')
    
    ps = numpy.zeros(len(samples), dtype=int)
    ps[:] = -2147483648 #set all to missing

    bases = vrecord.gt_bases.ravel()
    pos = vrecord.POS

    #walk across samples with phasing information
    for i in idx:

        #get alleles from GT field
        alleles = bases[i] #get the alleles
        if'|' in alleles:
            allele1, allele2 = alleles.split('|')
        else:
            allele1, allele2 = alleles.split('/')

        allele1_idx, allele2_idx, phase = genotypes[i]
        phase = False #default to unphased


        if allele1_idx == allele2_idx:
            stats['phase_not_transferred_for_homozygous_variant'] = stats.get('phase_not_transferred_for_homozygous_variant', 0) + 1
            genotypes[i] = [allele1_idx, allele2_idx, phase]
            continue

        #parse psw/gtw
        ps_id = psw[i]
        sample = samples[i]
        pos3, ref_len, phased_alleles, delimiter = parse_gtw(gtw[i])
        pallele1,pallele2 = phased_alleles

        vcf_has_star = '*' in alleles
        phased_has_star = any(['*' in x for x in phased_alleles])

        orientation = 0
        ref_len_match = ref_len == len(vrecord.REF)

        #start matching alleles
        phased_pos, pos_check = complete_pos(pos3, pos)


        #option 1: position matches, alleles match
        if pallele1 == allele1 and pallele2 == allele2 and pos_check:
            orientation = 1
            phase=True
            stats['phase_transferred'] = stats.get('phase_transferred', 0) + 1
        #option 2: position matches, alleles match in reverse
        elif pallele1 == allele2 and pallele2 == allele1 and pos_check:
            orientation = 2
            phase=True
            stats['phase_transferred']  = stats.get('phase_transferred', 0) + 1
        #option 3: position matches, possible suffix was added to phased alleles
        elif pos_check and not ref_len_match: #if the reference length is not the same, but position is, we need to drop the common suffix
            if not (vcf_has_star or phased_has_star):
                orientation, phase = check_common_suffix(allele1, allele2, pallele1, pallele2, ref_len, len(vrecord.REF), stats)
        #option 4: variant in VCF covers deletion in phasing
        elif phased_pos <= pos and vcf_has_star and not phased_has_star: #possible deletion in sample, covered by other variant in VCF
            orientation, phase = check_vcf_deletion(vrecord.REF, allele1, allele2, pallele1, pallele2, ref_len, pos, phased_pos, stats)
        elif vcf_has_star and phased_has_star: #both have deletions, check if they overlap
            orientation, phase = check_double_deletion(vrecord.REF, allele1, allele2, pallele1, pallele2, ref_len, pos, phased_pos, stats)
        else:
            pass

        #if the orientation was reversed, fix it here.
        if orientation == 2:
            allele1_idx, allele2_idx = allele2_idx, allele1_idx


        #if no phasing was transferred, let the user know.
        if not phase:
            warning = True
            if pos_check and ref_len_match:
                stats['phase_not_transferred_allele_mismatch']  = stats.get('phase_not_transferred_allele_mismatch', 0) + 1
                warning=False
            elif pos_check:
                stats['phase_not_transferred_with_correct_pos']  = stats.get('phase_not_transferred_with_correct_pos', 0) + 1
            else:
                stats['phase_not_transferred_with_incorrect_pos']  = stats.get('phase_not_transferred_with_incorrect_pos', 0) + 1

            last_failed[sample] = pos

            if warning:
                print(f"Warning: phasing information for sample {samples[i]} at position {pos} does not match alleles in GT field. Setting GT to unphased.")
                print(f"GTW: {gtw[i]}, {vrecord.REF}:{vrecord.ALT} alleles: {alleles}, phased_alleles: {phased_alleles}, pos_check: {pos_check}, {pos3}, ref_len_match: {ref_len_match}, {ref_len}, {len(vrecord.REF)}")
            

            
        #write out the genotypes
        genotypes[i] = [allele1_idx, allele2_idx, phase]
        if phase:
            ps[i] = ps_id

        #now match the alleles.

        #2. store in dict with key = reduced ref,alt comb and value = (ref,alt, npref, nsuf) tuple
        #3. drop common prefix and suffix from alleles
        #4. match alleles to dict in both directions. 
        #5. If matching, check npref and nsuf to see if this also matches. Also check if the position
        #    matches after correcting for the prefix.
        #6. if match found, set GT phase info and PS info. Otherwise, set GT to unphased and do not set PS info, and print warning.

    vrecord.set_format('PS', ps.reshape(psw.shape))
    vrecord.genotypes = genotypes

    
    return vrecord



def main():
    parser = argparse.ArgumentParser(description='Ligate and apply phasing information from merge_phasing.py after joint genotyping.')

    # Adding three required positional arguments
    parser.add_argument('--vcf', help='Input VCF')
    parser.add_argument('--output_vcf', help='Output VCF')
    parser.add_argument('--output_stats' , help='Output stats')

    # Adding an optional quiet flag
    parser.add_argument('-q', '--quiet', action='store_true', help='Run in quiet mode')

    args = parser.parse_args()

    # Open the input VCF files
    vcf = cyvcf2.VCF(args.vcf)
    samples = vcf.samples

    # Open the output gVCF file for writing
    out = cyvcf2.Writer(args.output_vcf, vcf)

    for file in [vcf,out]:
        file.add_format_to_header({'ID': 'PS', 'Number': '1', 'Type': 'Integer', 'Description': 'Phase block ID for phase info in GT.'})

    stats = {}
    write_stack = []
    
    last_failed = {s:0 for s in samples}
    sys.stderr.write('Processing variants...\n\n')
    try:
        vcf_record = next(vcf)
    except StopIteration:
        vcf_record = None
    while vcf_record is not None:
        stats['nvariant'] = stats.get('nvariant', 0) + 1
        sys.stderr.write(f"\r {stats.get('nvariant_with_gtw',0)}:{stats.get('nvariant',0)}")
        sys.stderr.flush()
        vcf_record = process_record(vcf_record, samples, stats, write_stack, last_failed)
        write_stack.append(vcf_record)

        if len(write_stack) > 200:
            for w in write_stack[:100]:
                out.write_record(w)
            write_stack = write_stack[100:]
        try:
            vcf_record = next(vcf)
        except StopIteration:
            vcf_record = None

    for w in write_stack:
        out.write_record(w)
    sys.stderr.write('Done processing ' + str(stats['nvariant']) + ' variants\n') 

    stats['nsamples'] = len(samples)

    # Close the output gVCF file
    out.close()

    with open(args.output_stats, 'w') as f:
        for k,v in stats.items():
            f.write(k + '\t' + str(v) + '\n')

    print(stats)

if __name__ == '__main__':
    main()
