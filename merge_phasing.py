import cyvcf2
import numpy
import sys


class PhaseBlock:
    def __init__(self, phase_id, source):
        self.variants = []
        self.source = source
        self.id = phase_id
        self.orig_id = phase_id
        self.correct_orientation = True
        self.parent=None
        self.children_blocks = set()
        self.invalid = False

    def addVariant(self, pos):
        self.variants.append(pos)

    def getAncestorBlock(self):
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

        
    def setParent(self, phase_block, same_orientation):
        
        if same_orientation is None:
            return False    
        ancestor_phase_block = phase_block.getAncestorBlock()
        if self is ancestor_phase_block:
            return False

        if ancestor_phase_block is not phase_block:        
            print(f'INFO: Connecting phase blocks {self} to {ancestor_phase_block} through {phase_block}')
        else:
            print(f'INFO: Connecting phase blocks {self} to {ancestor_phase_block}')            
        if isinstance(self.parent, PhaseBlock):
            if self.parent.getAncestorBlock() is ancestor_phase_block:
                if same_orientation is True: 
                    return True
                else:
                    print(f'ERROR: Phase block orientation mismatch, disconnecting phase blocks {self} and {ancestor_phase_block}')
                    stats['nphased_orientation_mismatch'] = stats.get('nphased_orientation_mismatch', 0) + 1
                    self.parent.removeChild(self)
                    self.parent = False
                    self.correct_orientation = True
                    self.id = self.orig_id
                    self.invalid = True
                    return False
            else:
                print(f'ERROR: Phase block ancestry mismatch, disconnecting phase blocks {self} and {ancestor_phase_block}')
                stats['nphased_ancestry_mismatch'] = stats.get('nphased_ancestry_mismatch', 0) + 1
                self.parent.removeChild(self)
                self.parent = False
                self.correct_orientation = True
                self.id = self.orig_id                                    
                self.invalid = True
                return False                
        else: 
            self.parent = phase_block
            self.parent.addChild(self)
            self.correct_orientation = same_orientation
            self.id = ancestor_phase_block.id
            return True
        
def process_record(grecord, vrecord, gphase_blocks, wphase_blocks):
    gphase_block = None
    if 'PID' in grecord.FORMAT:
        gatk_phase_block = int(grecord.format('PID')[0].split('_')[0])
        gphase_block = gphase_blocks[gatk_phase_block]
        gvariant = grecord.format('PGT').ravel()[0]
        gatk_allele1, gatk_allele2 =  [int(e) for e in gvariant.split('|')]
        gatk_allele1 = grecord.REF if gatk_allele1 == 0 else grecord.ALT[gatk_allele1-1]
        gatk_allele2 = grecord.REF if gatk_allele2 == 0 else grecord.ALT[gatk_allele2-1]
        gvariant = f'{gatk_allele1}|{gatk_allele2}'
        if gphase_block.invalid:
            gphase_block = None

    wphase_block = None        
    if vrecord is not None:
        whatshap_phase_block = int(vrecord.format('PS')[0])
        wphase_block = wphase_blocks[whatshap_phase_block]
        wvariant = vrecord.gt_bases.ravel()[0]
        if wphase_block.invalid:
            wphase_block = None

    
    if gphase_block is None and wphase_block is None:
        return grecord
    elif gphase_block is None:
        phase_block = wphase_block
        variant = wvariant
        vpos = vrecord.POS
        vreflen = len(vrecord.REF)
        delim = '|'
    elif wphase_block is None:
        phase_block = gphase_block
        variant = gvariant
        vpos = grecord.POS
        delim = '!'
        vreflen = len(grecord.REF)
    else:
        #if gphase_block.getAncestorBlock() is not  wphase_block.getAncestorBlock():
        #    stats['whatshap_gatk_ancestor_mismatch'] = stats.get('whatshap_gatk_ancestor_mismatch', 0) + 1
        #prefer whatshap
        phase_block = wphase_block
        variant = wvariant
        vpos = vrecord.POS 
        delim = '|'
        vreflen = len(vrecord.REF) if wphase_block.id < gphase_block.id else len(grecord.REF)

    block_id = phase_block.id
    v_allele1, v_allele2 = variant.split('|')
    if not phase_block.correct_orientation:
        v_allele1, v_allele2 = v_allele2, v_allele1
    vpos_str = str(vpos)[-3:]
    if vreflen == 1:
        vreflen = ''
    else:
        vreflen = f'+{vreflen}'
    result = f'{vpos_str}{vreflen}{v_allele1}{delim}{v_allele2}'
    
    grecord.set_format('GTW',numpy.array([result],dtype='S'))
    grecord.set_format('PSW', numpy.array([block_id],dtype='i'))                    
    return grecord
        
        
# Open the gVCF file
gvcf = cyvcf2.VCF(sys.argv[1])

# Open the phased VCF file
vcf = cyvcf2.VCF(sys.argv[2])

# Open the output gVCF file for writing
out = cyvcf2.Writer(sys.argv[3], gvcf)


# Add the PGT format field to the header of the output gVCF file

for file in [gvcf,out]:
    file.add_format_to_header({'ID': 'GTW', 'Number': '1', 'Type': 'String', 'Description': 'Alleles phased by a combination of Whatshap (separator |) and GATK physical phasing (separator !). Last 3 digits of position in GVCF, followed by + and the length of the reference allele if it is longer than 1 base, followed by the alleles separated by |.'})
    file.add_format_to_header({'ID': 'PSW', 'Number': '1', 'Type': 'Integer', 'Description': 'Phase block ID for phase info in GTW.'})

# Get the next record from each file
gvcf_record = next(gvcf)
vcf_record = next(vcf)
while 'PS' not in vcf_record.FORMAT:
    vcf_record = next(vcf)
write_stack = []
gphase_blocks = {}
wphase_blocks = {}

stats = {'nphased_variants_wh': 1}
# Loop over the variants in both files
while gvcf_record is not None and vcf_record is not None:
    # If the gVCF record is before the phased VCF record, get the next gVCF record
    if gvcf_record.CHROM < vcf_record.CHROM or (gvcf_record.CHROM == vcf_record.CHROM and gvcf_record.POS < vcf_record.POS):
        if 'PID' in gvcf_record.FORMAT:
            gatk_phase_block = int(gvcf_record.format('PID')[0].split('_')[0])
            if not gatk_phase_block in gphase_blocks:
                gphase_blocks[gatk_phase_block] = PhaseBlock(gatk_phase_block, 'GATK')
            gphase_blocks[gatk_phase_block].addVariant(gvcf_record.POS)
                     
        write_stack.append((gvcf_record, None))
        try:
            gvcf_record = next(gvcf)
        except StopIteration:
            gvcf_record = None
            
    # If the phased VCF record is before the gVCF record, get the next phased VCF record
    elif vcf_record.CHROM < gvcf_record.CHROM or (vcf_record.CHROM == gvcf_record.CHROM and vcf_record.POS < gvcf_record.POS):
        try:
            vcf_record = next(vcf)
            while 'PS' not in vcf_record.FORMAT:
                vcf_record = next(vcf)
                
            stats['nphased_variants_wh'] = stats.get('nphased_variants_wh',0) + 1
        except StopIteration:
            vcf_record = None
    # If the gVCF and phased VCF records position match, update the gVCF record with the phased genotype information
    else:
        if 'PID' in gvcf_record.FORMAT:
            gatk_phase_block = int(gvcf_record.format('PID')[0].split('_')[0])
            if not gatk_phase_block in gphase_blocks:
                gphase_blocks[gatk_phase_block] = PhaseBlock(gatk_phase_block, 'GATK')
            gphase_blocks[gatk_phase_block].addVariant(gvcf_record.POS)
            
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
            wh_allele1, wh_allele2 = vcf_record.gt_bases.ravel()[0].split("|")
            wh_phase_block = vcf_record.format('PS').ravel()[0]
            
            if not wh_phase_block in wphase_blocks:
                wphase_blocks[wh_phase_block] = PhaseBlock(wh_phase_block, 'whatshap')
            wphase_blocks[wh_phase_block].addVariant(gvcf_record.POS)
            
            if not wphase_blocks[wh_phase_block].correct_orientation:
                wh_allele1, wh_allele2 = wh_allele2, wh_allele1
                
            #two phase blocks GATK and whatshap overap
            if 'PID' in gvcf_record.FORMAT:              
                #determine if gatk and whatshap phase blocks are aligned
                gatk_allele1, gatk_allele2 =  [int(e) for e in gvcf_record.format('PGT').ravel()[0].split('|')]
                gatk_allele1 = gref if gatk_allele1 == 0 else galt[gatk_allele1-1]
                gatk_allele2 = gref if gatk_allele2 == 0 else galt[gatk_allele2-1]
                if not gphase_blocks[gatk_phase_block].correct_orientation:
                    gatk_allele1, gatk_allele2 = gatk_allele2, gatk_allele1
                    
                if wh_allele1 == gatk_allele1 and wh_allele2 == gatk_allele2:
                    same_orientation = True
                elif wh_allele2 == gatk_allele1 and wh_allele1 == gatk_allele2:
                    same_orientation = False
                else:
                    same_orientation = None
                    print('whatshap and gatk alleles do not match', gvcf_record.POS, wh_allele1, wh_allele2, gatk_allele1, gatk_allele2)
                    stats['whatshap_gatk_allele_mismatch'] = stats.get('whatshap_gatk_allele_mismatch',0) + 1
                    
                                          
                if not same_orientation is None:
                    if gphase_blocks[gatk_phase_block].isDominant(wphase_blocks[wh_phase_block]):
                        wphase_blocks[wh_phase_block].setParent(gphase_blocks[gatk_phase_block], same_orientation)
                    else:
                        gphase_blocks[gatk_phase_block].setParent(wphase_blocks[wh_phase_block], same_orientation)
                    
            # Write the updated record to the output gVCF file
            write_stack.append((gvcf_record, vcf_record))
            stats['nphased_variants_wh_linked'] = stats.get('nphased_variants_wh_linked',0) + 1                          
        else:
            # Write the updated record to the output gVCF file
            write_stack.append((gvcf_record, None))
        # Get the next record from each file
        try:
            gvcf_record = next(gvcf)
        except StopIteration:
            gvcf_record = None
        

            
            
    # Write the records in the write stack to the output gVCF file if the write stack is large        
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
        if not gatk_phase_block in gphase_blocks:
            gphase_blocks[gatk_phase_block] = PhaseBlock(gatk_phase_block, 'GATK')
        gphase_blocks[gatk_phase_block].addVariant(gvcf_record.POS)
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
s = set()
for v in variants:
    if v.intersection(s):
        print(v.intersection(s))
    s.update(v)

nvariants = len(s)
stats['nvariants_phased'] = nvariants
blocklength = [max(v) - min(v) for v in variants]
stats['average_phase_block_length'] = numpy.mean(blocklength)
stats['median_phase_block_length'] = numpy.median(blocklength)
stats['max_phase_block_length'] = max(blocklength)
stats['sum_phase_block_length'] = sum(blocklength)


# Close the output gVCF file
out.close()

with open(sys.argv[4], 'w') as f:
    for k,v in stats.items():
        f.write(k + '\t' + str(v) + '\n')

print(stats)        