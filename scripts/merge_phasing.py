import argparse

import cyvcf2
import numpy

stats = {}

class PhaseBlock:
    """Represents a phase block. A phase block is a set of variants that are phased together."""
    def __init__(self, phase_id, source):
        self.variants = []
        self.source = source
        self.id = phase_id #block id, format is chrom:start position of the block. Can be changed if the block is connected to another block to the ancestor id of that block
        self.orig_id = phase_id
        self.correct_orientation = True #orrientation wrt ancestor block. If False, the alleles should be swapped
        self.parent=None #block with a lower id (i.e. started earlier on the chromosome) which overlaps with this block
        self.children_blocks = set()
        self.invalid = False #if the block is invalid, there was a phasing mismatch, and the block should not be used. This is set to True if the block is connected to conflicting ancestors.

    def addVariant(self, pos):
        if self.variants:
            assert pos >= self.variants[-1], 'Adding variant with earlier position should not be possible'
            
            if (pos - self.variants[-1]) > 50000: #check if the block makes too big of a jump. Happens rarely, should be filtered. 
                print(f"Jump attempted for {self} by adding variant {pos} to {self.variants}")
                return False
        self.variants.append(pos)
        return True

    def getAncestorBlock(self):
        #look up the parent of the parent until you find a block without a parent
        if isinstance(self.parent, PhaseBlock):
            return self.parent.getAncestorBlock()
        else:
            return self

    def isDominant(self, other):
        if self.id < other.id:
            return True
        else:
            return False

    def removeChild(self, phase_block):
        self.children_blocks.remove(phase_block)

    def addChild(self, phase_block):
        self.children_blocks.add(phase_block)

    def getAllAncestorBlocks(self):
        if isinstance(self.parent, PhaseBlock):
            return [self.parent] + self.parent.getAllAncestorBlocks()
        else:
            return []

    def __hash__(self):
        return hash((self.orig_id, self.source))

    def __eq__(self, other):
        return self.orig_id == other.orig_id and self.source == other.source

    def __lt__(self, other):
        return self.orig_id < other.orig_id or (self.orig_id == other.orig_id and self.source < other.source)

    def fullPrint(self, prefix = ''):
        res = f'{prefix}PB-{self.source} {self.orig_id}: {self.variants}'
        for child in self.children_blocks:
            cp = child.fullPrint(prefix + "\t - ")
            res += f'\n{cp}'
        return res

    def getAllVariants(self):
        res = self.variants
        for child in self.children_blocks:
            res += child.getAllVariants()
        return set(res)

    def __repr__(self):
        ancestors = '>'.join([f'{e.source}:{e.orig_id}' for e in self.getAllAncestorBlocks()])

        res = f'PB-{self.source} {self.orig_id}'

        if ancestors:
            res += f'>{ancestors}'
        res += f': {self.variants}'
        return res


    def setParent(self, phase_block, same_orientation, quiet):
        """Set the parent of this phase block to the given phase block. 
        If the parent is already set, check if the parent is the same as the given phase block. 
        If not, disconnect the parent and set the parent to the given phase block.
        
        Same_orientation is True if the orientation of this block is the same as the orientation of the parent block. Should be w.r.t. to the corrected orientation of the parent block and this block.
        """
        if same_orientation is None:
            return False
        ancestor_phase_block = phase_block.getAncestorBlock()
        if self is ancestor_phase_block:
            return False
        if not quiet:
            if ancestor_phase_block is not phase_block:
                print(f'INFO: Connecting phase blocks {self} to {ancestor_phase_block} through {phase_block}')
            else:
                print(f'INFO: Connecting phase blocks {self} to {ancestor_phase_block}')
        if isinstance(self.parent, PhaseBlock): #if this block already has a parent
            if self.parent.getAncestorBlock() is ancestor_phase_block: #if the parent of this block is the same as the parent of the block we want to connect to
                
                if same_orientation is True: 
                    #is this correct? or should we check self.correct_orientation == same_orientation? or even combine the parent correct_orientations?
                    # --> this is correct. correct_orrientation is w.r.t. to the ancestor block, so if the ancestor blocks are the same, the orientations should be the same.
                    #  --> also, correct_orientation for this block is already taken into account when calculating same_orientation. Correct_orientation is wrt ancestor,
                    #       and correct_orrientation of the parent also wr.t. to the same ancestor, so same_orientation needs to be True. 
                    
                    return True #nothing to do
                
                    #if this is not the case, this is an error.
                else:
                    #
                    print(f'ERROR: Phase block orientation mismatch with {phase_block}, disconnecting phase blocks {self} and {ancestor_phase_block}')
                    stats['nphased_orientation_mismatch'] = stats.get('nphased_orientation_mismatch', 0) + 1
                    self.parent.removeChild(self)
                    self.parent = False
                    self.correct_orientation = True
                    self.id = self.orig_id
                    self.invalid = True
                    return False
            else:
                #we have to set the parent, but the parent is not compatible with the parent we already have.
                #that means that the (ancestor of the) parent we already have is different from the (ancestor of the) parent we want to set
                #this means that the parents are not connected. This is weird, as in read-based phasing, if they overlap, they should be connected
                print(f'ERROR: Phase block ancestry mismatch over {phase_block}, disconnecting phase blocks {self} and {ancestor_phase_block}')
                stats['nphased_ancestry_mismatch'] = stats.get('nphased_ancestry_mismatch', 0) + 1
                self.parent.removeChild(self)
                self.parent = False
                self.correct_orientation = True
                self.id = self.orig_id
                self.invalid = True
                return False
        else: # no parent yet, so set it. This is the most common case
            self.parent = phase_block
            self.parent.addChild(self)
            #set the orientation of this block to the same as the parent. 
            # Should we not check the orientation of the parent? --> no, same_orientation is wrt. to 
            # corrected orientation of the parent, see below
            # (and this block shold have correct_orientation = True as is did not have a parent yet)
            self.correct_orientation = same_orientation 
            self.id = ancestor_phase_block.id
            return True

def process_record(grecord, vrecord, gphase_blocks, wphase_blocks):
    """Process a gVCF record and a phased VCF record (if available at this position) and return the updated gVCF record"""
    gphase_block = None
    if 'PID' in grecord.FORMAT: #if gvcf record is phased
        gatk_phase_block = int(grecord.format('PID')[0].split('_')[0])
        gatkkey = f'{grecord.CHROM}:{gatk_phase_block}'
        gphase_block = gphase_blocks[gatkkey]
        gvariant = grecord.format('PGT').ravel()[0]
        gatk_allele1, gatk_allele2 =  (int(e) for e in gvariant.split('|')) #get the alleles
        gatk_allele1 = grecord.REF if gatk_allele1 == 0 else grecord.ALT[gatk_allele1-1]
        gatk_allele2 = grecord.REF if gatk_allele2 == 0 else grecord.ALT[gatk_allele2-1]
        gvariant = f'{gatk_allele1}|{gatk_allele2}'
        if gphase_block.invalid: #if the phase block is invalid, we should not use it
            gphase_block = None 

    # Get the phased genotype information from the VCF record
    wphase_block = None
    if vrecord is not None:
        whatshap_phase_block = int(vrecord.format('PS')[0])
        whkey = f'{vrecord.CHROM}:{whatshap_phase_block}'
        wphase_block = wphase_blocks[whkey]
        wvariant = vrecord.gt_bases.ravel()[0]
        if wphase_block.invalid:
            wphase_block = None


    if gphase_block is None and wphase_block is None: #if neither gvcf nor vcf record is phased
        return grecord
    elif gphase_block is None: #if only the vcf record is phased
        phase_block = wphase_block #use the whatshap phase block (if it is not invalid)
        variant = wvariant
        vpos = vrecord.POS
        vreflen = len(vrecord.REF)
        delim = '|'
    elif wphase_block is None: #if only the gvcf record is phased
        phase_block = gphase_block #use the gatk phase block (if it is not invalid)
        variant = gvariant
        vpos = grecord.POS
        delim = '!'
        vreflen = len(grecord.REF)
    else: #if both gvcf and vcf record are phased
        #if gphase_block.getAncestorBlock() is not  wphase_block.getAncestorBlock():
        #    stats['whatshap_gatk_ancestor_mismatch'] = stats.get('whatshap_gatk_ancestor_mismatch', 0) + 1
        #prefer whatshap
        phase_block = wphase_block #prefer whatshap
        variant = wvariant
        vpos = vrecord.POS
        delim = '|'
        vreflen = len(vrecord.REF) 

    block_id = int(str(phase_block.id).split(':')[1]) #output phase block id
    v_allele1, v_allele2 = variant.split('|')  #get the alleles
    if not phase_block.correct_orientation: #if the phase block is in the wrong orientation, swap the alleles
        v_allele1, v_allele2 = v_allele2, v_allele1
    vpos_str = str(vpos)[-3:] #last 3 digits of position
    if vreflen == 1: #if the reference allele is 1 base, do not output the length
        vreflen = ''
    else: #if the reference allele is longer than 1 base, output the length
        vreflen = f'+{vreflen}' #why? can we not just derive this from the alleles --> to clarify which one is the ref allele. 
    result = f'{vpos_str}{vreflen}{v_allele1}{delim}{v_allele2}'

    grecord.set_format('GTW',numpy.array([result],dtype='S'))
    grecord.set_format('PSW', numpy.array([block_id],dtype='i'))
    return grecord



def main():
    parser = argparse.ArgumentParser(description='Merge WhatsHap phasing information into GATK GVCF')

    # Adding three required positional arguments
    parser.add_argument('gvcf', help='Input GVCF')
    parser.add_argument('vcf',  help='Phased VCF')
    parser.add_argument('output_gvcf', help='Output GVCF')
    parser.add_argument('output_stats' , help='Output stats')

    # Adding an optional quiet flag
    parser.add_argument('-q', '--quiet', action='store_true', help='Run in quiet mode')

    args = parser.parse_args()
    # Open the gVCF file
    gvcf = cyvcf2.VCF(args.gvcf)

    # Open the phased VCF file
    vcf = cyvcf2.VCF(args.vcf)

    # Open the output gVCF file for writing
    out = cyvcf2.Writer(args.output_gvcf, gvcf)


    # Add the PGT format field to the header of the output gVCF file

    for file in [gvcf,out]:
        file.add_format_to_header({'ID': 'GTW', 'Number': '1', 'Type': 'String', 'Description': 'Alleles phased by a combination of Whatshap (separator |) and GATK physical phasing (separator !). Last 3 digits of position in GVCF, followed by + and the length of the reference allele if it is longer than 1 base, followed by the alleles separated by |.'})
        file.add_format_to_header({'ID': 'PSW', 'Number': '1', 'Type': 'Integer', 'Description': 'Phase block ID for phase info in GTW.'})


    # Get the next record from each file
    try:
        gvcf_record = next(gvcf)
    except StopIteration:
        gvcf_record = None

    try:
        vcf_record = next(vcf)
        while 'PS' not in vcf_record.FORMAT:
            vcf_record = next(vcf)
        stats['nphased_variants_wh'] = 1
    except StopIteration:
        vcf_record = None
        stats['nphased_variants_wh'] = 0

    write_stack = []
    gphase_blocks = {}
    wphase_blocks = {}

    # Loop over the variants in both files
    while gvcf_record is not None and vcf_record is not None:
        # If the gVCF record is before the phased VCF record, get the next gVCF record
        if gvcf_record.CHROM < vcf_record.CHROM or (gvcf_record.CHROM == vcf_record.CHROM and gvcf_record.POS < vcf_record.POS):
            if 'PID' in gvcf_record.FORMAT: #if gvcf record is phased
                gatk_phase_block = int(gvcf_record.format('PID')[0].split('_')[0]) #get the phase block id
                gatkkey = f'{gvcf_record.CHROM}:{gatk_phase_block}'
                if gatkkey not in gphase_blocks: #is this a new phase block?
                    gphase_blocks[gatkkey] = PhaseBlock(gatkkey, 'GATK')
                gphase_blocks[gatkkey].addVariant(gvcf_record.POS) #add the variant to the phase block

            write_stack.append((gvcf_record, None)) 
            try:
                gvcf_record = next(gvcf)
            except StopIteration:
                gvcf_record = None

        # If the phased VCF record is before the gVCF record, get the next phased VCF record
        elif vcf_record.CHROM < gvcf_record.CHROM or (vcf_record.CHROM == gvcf_record.CHROM and vcf_record.POS < gvcf_record.POS):
            try:
                vcf_record = next(vcf)
                while 'PS' not in vcf_record.FORMAT: #skip unphased variants
                    vcf_record = next(vcf)

                stats['nphased_variants_wh'] = stats.get('nphased_variants_wh',0) + 1
            except StopIteration:
                vcf_record = None
        else: #if the records are at the same position
            if 'PID' in gvcf_record.FORMAT: #if gvcf record is phased
                gatk_phase_block = int(gvcf_record.format('PID')[0].split('_')[0]) #get the phase block id
                gatkkey = f'{gvcf_record.CHROM}:{gatk_phase_block}'
                if gatkkey not in gphase_blocks: #is this a new phase block?
                    gphase_blocks[gatkkey] = PhaseBlock(gatkkey, 'GATK') 
                gphase_blocks[gatkkey].addVariant(gvcf_record.POS) #add the variant to the phase block

            # Check if the reference and alternate alleles match
            # FIXME: this does not account for the case where the gVCF record has a longer reference allele than the phased VCF record
            # Of 5342 records, this matches 5258, so it is not a big problem
            # Two options to account for:
            # 1. position is the same, but reference allele is longer in gVCF record
            # 2. position is also different
            # lets first try to fix option 1.
            gref = gvcf_record.REF
            vref = vcf_record.REF
            galt = gvcf_record.ALT
            valt = vcf_record.ALT
            if len(gref) > len(vref) and gref.startswith(vref):
                add_to_alt = gref[len(vref):]
                valt = [alt + add_to_alt for alt in valt]
                vref = gref
            #this fixes all records in my testing (5342/5342 matched)

            if gref == vref and set(valt).issubset(set(galt)):
                

                # Get the phased genotype information from the VCF record
                assert vcf_record.gt_phases, 'PS field present without phasing information'
                wh_allele1, wh_allele2 = vcf_record.gt_bases.ravel()[0].split("|") #get the alleles
                wh_phase_block = vcf_record.format('PS').ravel()[0] #get the phase block id
                whkey = f'{vcf_record.CHROM}:{wh_phase_block}'

                if whkey not in wphase_blocks: #is this a new phase block?
                    wphase_blocks[whkey] = PhaseBlock(whkey, 'whatshap')
                if wphase_blocks[whkey].addVariant(gvcf_record.POS):
                    #re-orient alleles if necessary (occurs if the phase block is connected to a parent block with a different orientation)
                    if not wphase_blocks[whkey].correct_orientation:
                        wh_allele1, wh_allele2 = wh_allele2, wh_allele1

                    #two phase blocks GATK and whatshap overlap
                    if 'PID' in gvcf_record.FORMAT:
                        #determine if gatk and whatshap phase blocks are aligned
                        gatk_allele1, gatk_allele2 =  (int(e) for e in gvcf_record.format('PGT').ravel()[0].split('|')) #get the alleles
                        gatk_allele1 = gref if gatk_allele1 == 0 else galt[gatk_allele1-1] #from allele index to allele
                        gatk_allele2 = gref if gatk_allele2 == 0 else galt[gatk_allele2-1]
                        #re-orient alleles if necessary
                        if not gphase_blocks[gatkkey].correct_orientation:
                            gatk_allele1, gatk_allele2 = gatk_allele2, gatk_allele1

                        if wh_allele1 == gatk_allele1 and wh_allele2 == gatk_allele2:
                            same_orientation = True
                        elif wh_allele2 == gatk_allele1 and wh_allele1 == gatk_allele2:
                            same_orientation = False
                        else:
                            same_orientation = None
                            print('whatshap and gatk alleles do not match', gvcf_record.POS, wh_allele1, wh_allele2, gatk_allele1, gatk_allele2)
                            stats['whatshap_gatk_allele_mismatch'] = stats.get('whatshap_gatk_allele_mismatch',0) + 1


                        if same_orientation is not None:
                            if gphase_blocks[gatkkey].isDominant(wphase_blocks[whkey]): #if gatk phase block is dominant (i.e. has a lower id, meaning it started earlier on the chromosome)
                                wphase_blocks[whkey].setParent(gphase_blocks[gatkkey], same_orientation,args.quiet) #set the parent of the whatshap phase block to the gatk phase block
                            else:
                                gphase_blocks[gatkkey].setParent(wphase_blocks[whkey], same_orientation,args.quiet) #else set the parent of the gatk phase block to the whatshap phase block

                    # Write the updated record to the output gVCF file
                    write_stack.append((gvcf_record, vcf_record))
                    stats['nphased_variants_wh_linked'] = stats.get('nphased_variants_wh_linked',0) + 1
                else:
                    write_stack.append((gvcf_record, None))
            else:
                # Write the updated record to the output gVCF file
                write_stack.append((gvcf_record, None))
            # Get the next record from each file
            try:
                gvcf_record = next(gvcf)
            except StopIteration:
                gvcf_record = None




        # Write the records in the write stack to the output gVCF file if the write stack is large
        # We only do this after processing 25000 records, to make sure the phase blocks are not rearranged.
        if(len(write_stack) > 25000):
            to_write = write_stack[:-24000]
            write_stack = write_stack[-24000:]
            for grecord, vrecord in to_write:
                grecord = process_record(grecord, vrecord, gphase_blocks, wphase_blocks)
                out.write_record(grecord)


    # Write the remaining records to the output gVCF file
    while gvcf_record is not None:
        if 'PID' in gvcf_record.FORMAT:
            gatk_phase_block = int(gvcf_record.format('PID')[0].split('_')[0])
            gatkkey = f'{gvcf_record.CHROM}:{gatk_phase_block}'
            if gatkkey not in gphase_blocks:
                gphase_blocks[gatkkey] = PhaseBlock(gatkkey, 'GATK')
            gphase_blocks[gatkkey].addVariant(gvcf_record.POS)
        write_stack.append((gvcf_record,None))
        try:
            gvcf_record = next(gvcf)
        except StopIteration:
            gvcf_record = None

    for grecord, vrecord in write_stack:
        grecord = process_record(grecord, vrecord, gphase_blocks, wphase_blocks)
        out.write_record(grecord)

    stats['nblocks_total'] = len(gphase_blocks) + len(wphase_blocks)
    stats['gatk_blocks'] = len(gphase_blocks)
    stats['whatshap_blocks'] = len(wphase_blocks)

    #bit of printing
    blocks = list(gphase_blocks.items()) + list(wphase_blocks.items())
    blocks = [(a,b) for a,b in blocks if b.parent is None]
    blocks.sort()

    print('All full blocks:')
    for bid, block in blocks:
        print(block.fullPrint())



    #some stats
    stats['nblocks_final'] = len(blocks)
    stats['whatshap_blocks_final'] = len([b for b in blocks if b[1].source == 'whatshap'])
    stats['gatk_blocks_final'] = len([b for b in blocks if b[1].source == 'GATK'])
    variants = [b[1].getAllVariants() for b in blocks]

    #detect variants in multiple blocks
    s = set()
    for v in variants:
        if v.intersection(s):
            print(v.intersection(s))
        s.update(v)

    nvariants = len(s)
    stats['nvariants_phased'] = nvariants
    blocklength = [max(v) - min(v) for v in variants]
    idx = numpy.argmax(blocklength)
    print('BLOCKLENGTH MAX', idx, blocklength[idx], variants[idx])
    print(blocks[idx])
    stats['average_phase_block_length'] = numpy.mean(blocklength)
    stats['median_phase_block_length'] = numpy.median(blocklength)
    if blocklength:
        stats['max_phase_block_length'] = max(blocklength)
        stats['sum_phase_block_length'] = sum(blocklength)
        stats['average_nvariants_per_block'] = nvariants / len(blocklength)
        stats['median_nvariants_per_block'] = numpy.median([len(v) for v in variants])
        stats['max_nvariants_per_block'] = max([len(v) for v in variants])
    else:
        stats['max_phase_block_length'] = 0
        stats['sum_phase_block_length'] = 0
        stats['average_nvariants_per_block'] = 0
        stats['median_nvariants_per_block'] = 0
        stats['max_nvariants_per_block'] = 0



    # Close the output gVCF file
    out.close()

    with open(args.output_stats, 'w') as f:
        for k,v in stats.items():
            f.write(k + '\t' + str(v) + '\n')

    print(stats)

if __name__ == '__main__':
    main()
