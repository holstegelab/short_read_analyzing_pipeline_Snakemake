import csv
import os
import numpy
import yaml


CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

def write_tsv(filename, header,data):
    with open(filename, 'w') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(header)
        for row in data:
            w.writerow(row)

def read_artifacts(filename):
    results = {}
    with open(filename, 'r') as f:
        for row in f:
            if row.startswith('#') or row.strip() == "":
                continue
            elif row.startswith('SAMPLE'):
                continue
            else:
                row = row.split("\t")
                sample, library, ref, alt, qscore = row[:5]
                qscore = float(qscore)
                key = (ref, alt)
                if key in results:
                    results[key] = min(results[key], qscore)
                else:
                    results[key] = qscore

    return [results[('A','C')], results[('A','G')], results[('A','T')],
            results[('C','A')], results[('C','G')], results[('C','T')],
            results[('G','A')], results[('G','C')], results[('G','T')],
            results[('T','A')], results[('T','C')], results[('T','G')]]
         

def read_detail_artifacts(filename):
    results = {}
    with open(filename, 'r') as f:
        for row in f:
            if row.startswith('#') or row.strip() == "":
                continue
            elif row.startswith('SAMPLE'):
                continue
            else:
                row = row.split("\t")
                if len(row ) == 11: #pre_adapter
                    sample, library, ref, alt, context, pro_ref_bases, pro_alt_bases, con_ref_bases, con_alt_bases, error_rate, q_score = row
                    value = (int(pro_ref_bases), int(pro_alt_bases), int(con_ref_bases), int(con_alt_bases))
                elif len(row) == 13: #bait bias
                    sample, library, ref, alt, context, fwd_ref_bases, fwd_alt_bases, rev_ref_bases, rev_alt_bases, fwd_error_rate, rev_error_rate, error_rate, q_score = row
                    value = (int(fwd_ref_bases), int(fwd_alt_bases), int(rev_ref_bases), int(rev_alt_bases))
                else:
                    raise RuntimeError('Unknown format')
                key = (ref.lower(), alt.lower(), (context[0] + context[2]).lower())
                
                if not key in results:
                    results[key] = value
                else:
                    value_old = results[key]
                    value = tuple([v1 + v2 for v1,v2 in zip(value, value_old)])
                    results[key] = value

    return results
   


def read_stats(filename):
    results = {}
    with open(filename,'r') as f:
        for row in f:
            if row.startswith('#') or row.startswith('CHK'):
                continue
            elif row.startswith('SN'):
                row = row.strip().split('\t')
                results[row[1].rstrip(':')] = float(row[2])
            else:
                break
    return results

def read_stats_full(filename):
    results = {}
    cycle_count = []
    ffq = []
    lfq = []
    mismatches = []
    n_count = []

    gcf = []
    gcl = []
    gcc = []
    insertsize = []
    rl = []
    indel = []
    ic = []
    cov = []
    gcd = []
    print(filename)
    with open(filename,'r') as f:
        for row in f:
            print(row)
            if row.startswith('#') or row.startswith('CHK'):
                continue
            elif row.startswith('SN'):
                row = row.strip().split('\t')
                results[row[1].rstrip(':')] = float(row[2])
            elif row.startswith('FFQ'):
                line = row[3:].strip()
                w = tuple([int(e) for e in line.split('\t')])
                total = float(sum(w[1:]))
                if total > 0:
                    ffq.append(sum([pos * float(e) for pos,e in enumerate(w[1:])]) / total)
                    cycle_count.append(total)
            elif row.startswith('LFQ'):
                line = row[3:].strip()
                w = tuple([int(e) for e in line.split('\t')])
                total = float(sum(w[1:]))
                if total > 0:
                    lfq.append(sum([pos * e for pos,e in enumerate(w[1:])]) / total)
            elif row.startswith('MPC'):
                line = row[3:].strip()
                w = tuple([int(e) for e in line.split('\t')])
                n_count.append(w[1])
                mismatches.append(sum(w[2:]))
            elif row.startswith('GCF'):
                line = row[3:].strip()
                gcpercentage, count = line.split('\t')
                gcf.append(tuple([float(gcpercentage), int(count)]))
            elif row.startswith('GCL'):
                line = row[3:].strip()
                gcpercentage, count = line.split('\t')
                gcl.append(tuple([float(gcpercentage), int(count)]))
            elif row.startswith('GCC'):
                line = row[3:].strip()
                line = line.split('\t')
                gcc.append(tuple([int(line[0])] + [float(e) for e in line[1:]]))
            elif row.startswith('IS'):
                line = row[2:].strip()
                insertsize.append(tuple([int(e) for e in line.split('\t')]))
            elif row.startswith('RL'):
                line = row[2:].strip()
                rl.append(tuple([int(e) for e in line.split('\t')]))
            elif row.startswith('ID'):
                line = row[2:].strip()
                indel.append(tuple([int(e) for e in line.split('\t')]))
            elif row.startswith('IC'):
                line = row[2:].strip()
                ic.append(tuple([int(e) for e in line.split('\t')]))
            elif row.startswith('COV'):
                line = row[3:].strip()
                line = line.split('\t')
                unbounded =  '<' in line[0]
                cov.append(tuple([int(unbounded)]  + [int(e) for e in line[1:]]))
            elif row.startswith('GCD'):
                line = row[3:].strip()
                gcd.append(tuple([float(e) for e in line.split('\t')]))
            else:
                continue


        results['qual_by_cycle'] = '|'.join(['%.1f,%.1f,%.4f,%.4f' % (a,b,c/float(e), d/float(e)) for a,b,c,d,e in zip(ffq, lfq, n_count, mismatches,cycle_count)])
        results['gc_read1'] = sum([a * b for a,b in gcf]) / max(float(sum([b for a,b in gcf])),1e-6)
        results['gc_read2'] = sum([a * b for a,b in gcl]) / max(float(sum([b for a,b in gcl])),1e-6)
        results['base_fractions'] = '|'.join(['%.2f,%.2f,%.2f,%.2f' % (a,c,g,t) for cycle,a,c,g,t,n,o in gcc])
        total_pairs = max(float(sum([b for a,b,c,d,e in insertsize])),1e-6)
        print("IS", insertsize)
        frow = insertsize[0]
        lrow = insertsize[-1]
        results['insert_maxsize'] = lrow[0]
        results['pair_inward_maxsize'] = (lrow[2] / total_pairs)
        results['pair_outward_maxsize'] = (lrow[3] / total_pairs)
        results['pair_other_maxsize'] = (lrow[4] / total_pairs)
        results['pair_outward_size0'] =(frow[3] / total_pairs)
        results['pair_other_size0'] = (frow[4] / total_pairs)
        nr_insertions = sum([i for length,i,d in indel])
        nr_deletions = sum([d for length,i,d in indel])
        results['insertions_frac'] = nr_insertions / float(results['bases mapped'])
        results['deletions_frac'] = nr_deletions / float(results['bases mapped'])
        results['insertions_avg'] = sum([length * i for length,i,d in indel]) / float(nr_insertions)
        results['deletions_avg'] = sum([length * d for length,i,d in indel]) / float(nr_deletions)
        region = sum([c for a,b,c in cov])
        results['coverage1'] = cov[0][2] / float(region)
        results['coverage1000'] = cov[-1][2] / float(region)
        g35 =  [d50 for g,s,d10,d25,d50,d75,d90 in gcd if g == 35.0]
        g50 = [d50 for g,s,d10,d25,d50,d75,d90 in gcd if g == 50.0]
        if len(g35) > 0:
            g35 = g35[-1]
        else:
            g35 = numpy.nan

        if len(g50) > 0:
            g50 = g50[-1]
        else:
            g50 = numpy.nan

        results['gc_depth35'] = g35
        results['gc_depth50'] = g50

    return results



def read_contam(filename):
    with open(filename,'r') as f:
        c = csv.reader(f, delimiter='\t')
        c = iter(c)
        headers = c.__next__()
        data = c.__next__()
        freemix = data[6]
    return freemix

def combine_sex_stats(samples, result_files, sex_reported):
    result = []
    header = ['sample', 'reported', 'sex', 'yyratio', 'auto', 'chrx', 'chry', 'intercept', 'slope', 'rvalue', 'xratio', 'yratio', 'dist_auto', 'dist_chrx', 'dist_chrx_ratio', 'dist_chry', 'dist_chry_ratio']
    for sample, result_file, reported in zip(samples,result_files, sex_reported):
            with open(result_file,'r') as f:
                x = yaml.load(f, Loader=yaml.FullLoader)
                for k in x.keys():
                    if k.startswith('dist'):
                        x[k] = ';'.join(['%.2f' % e for e in x[k]])
                v = tuple([sample, reported] + [x[h] for h in header[2:]])                
                result.append(v)
    return (header,result)


def read_hsstats(filename):
    header = []
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
    while lines:
        row = lines.pop(0).strip()
        if 'BAIT_SET' in row:
            header = tuple([e.lower() for e in row.split('\t')])[:-2]
        elif header:
            data = tuple(row.split("\t"))
            break

    data= dict(zip(header,data))
    header2 = []
    data2 = []
    while lines:
       row = lines.pop(0).strip()
       if 'coverage' in row:
           header2 = row.split('\t')
       elif header2:
           if row:
               data2.append(tuple(row.split('\t')))

    
    return (data, data2)

def read_adapter_removal_stats(filename):
    result = {}
    with open(filename,'r') as f:
        rows = f.readlines()
        while rows:
            row = rows.pop(0)
            if 'Trimming statistics' in row:
                break   
        while rows:
            row = rows.pop(0)
            if ':' in row:   
                a,b = row.split(':')
                a = a.strip()
                if a.startswith('Number of'):
                    a = a[10:]
                elif a == 'Total number of read pairs':
                    a = 'ar_read_pairs'
                elif a == 'Average length of retained reads':
                    a = 'ar_average_read_length_retained'
                a = a.replace(' ', '_')
                b = float(b.strip())
                result[a] = b
            else:
                break 
    
    result['ar_aligned_fraction'] = result['well_aligned_read_pairs'] / result['ar_read_pairs']
    result['ar_adapter_fraction'] = result['reads_with_adapters[1]']  / (2.0 * result['ar_read_pairs'])
    result['ar_retained_fraction'] = result['retained_reads'] / (2.0 * result['ar_read_pairs'])
    return result 
    

def read_adapter_identify_stats(filename):
    with open(filename,'r') as f:
        rows = f.readlines()
        overlapping = int(rows[1][9:].strip().split(' ')[0])
        contained = int(rows[2][12:].strip().split(' ')[0])
        cons = [e[15:].strip() for e in rows if 'Consensus' in e]
    

    return {'ai_reads_aligned': overlapping, 'ai_n_adapter': contained, 'ai_adapter1':cons[0], 'ai_adapter2':cons[1]}


def read_merge_stats(filename):
    result = {}
    with open(filename,'r') as f:
        rows = f.readlines()
    while rows:
        row = rows.pop(0).strip()
        if not row:
            break 
        key,value = row.split('\t')
        result['merge_' + key] = float(value)
    result['merge_readded_fragments_ratio'] = result['merge_readded_fragments'] / result['merge_fragments']
    result['merge_supplementary_alignments_ratio'] = result['merge_supplementary_alignments'] / result['merge_alignments']
    result['merge_restored_read_ratio'] = (result['merge_restored_read1s'] + result['merge_restored_read2s']) / (2.0 * result['merge_fragments'])
    result['merge_restored_bp_ratio'] = (result['merge_restored_bp_read1'] + result['merge_restored_bp_read2']) / result['merge_total_bp']
    return result

def read_dechimer_stats(filename):
    result = {}
    with open(filename,'r') as f:
        for row in f.readlines():
            key,value = row.split('\t')
            result['dechimer_' + key] = float(value)
    result['dechimer_fragment_modified_ratio'] = result.get('dechimer_fragment_modified',0) / result.get('dechimer_fragment_counter',1.0)
    result['dechimer_clip_ratio'] = (result.get('dechimer_read1_dechimer_clip',0) + result.get('dechimer_read2_dechimer_clip',0)) / (2.0 * result.get('dechimer_fragment_counter',1.0))
    result['dechimer_loose_end_clip_ratio'] = (result.get('dechimer_read1_loose_end_clip',0) + result.get('dechimer_read2_loose_end_clip',0)) / (2.0 * result.get('dechimer_fragment_counter',1.0))
    result['dechimer_all_sup_discarded_ratio'] = (result.get('dechimer_read1_all_sup_discarded_in_pruning',0) + result.get('dechimer_read2_all_sup_discarded_in_pruning',0)) / \
                                                float(result.get('dechimer_read1_has_supplementary_alignments',1.0) + result.get('dechimer_read2_has_supplementary_alignments',1.0))
    result['dechimer_min_alignment_length_unmap_ratio'] = (result.get('dechimer_read1_min_align_length_unmap',0) + result.get('dechimer_read2_min_align_length_unmap',0)) / \
                                        (2.0 * result.get('dechimer_fragment_counter',1.0))
        
    return result
def valid_name(name):
    name = name.lower()
    newname = []
    for char in name:
        if(char.isalnum() or char == "_"):
            newname.append(char)
        else:
            newname.append("_")
    newname = "".join(newname)
    newname = newname.strip("_")
    if(not newname):
        return "funknown"
    if(newname[0].isdigit()):
        newname = "_" + newname
    while '__' in newname:
        newname = newname.replace('__', '_')
        
    return str(newname)

def read_dragmap_stats(filename):
    result_values = {}
    result_ratio = {}
    with open(filename, 'r') as f:
        rows = f.readlines()
    
    mappingrows = [e.split(',')[2:] for e in rows if 'MAPPING' in e]
    for row in mappingrows:
        key = 'dm_' + valid_name(row[0].strip())
        value = float(row[1].strip())
        if len(row) > 2:
            ratio = float(row[2].strip())
            result_ratio[key] = ratio
        result_values[key] = value

    result_ratio['dm_reads_with_indel'] = (result_ratio['dm_reads_with_indel_r1'] + result_ratio['dm_reads_with_indel_r2']) / 2.0        
    result_ratio['dm_soft_clipped_bases'] = (result_ratio['dm_soft_clipped_bases_r1'] + result_ratio['dm_soft_clipped_bases_r2']) / 2.0
    result_ratio['dm_mismatched_bases'] = (result_ratio['dm_mismatched_bases_r1'] + result_ratio['dm_mismatched_bases_r2']) / 2.0
    result_ratio['dm_mismatched_bases_excl_indels'] = (result_ratio['dm_mismatched_bases_r1_excl_indels'] + result_ratio['dm_mismatched_bases_r2_excl_indels']) / 2.0
    result_ratio['dm_q30_bases_diff'] = (result_ratio['dm_q30_bases_r1'] - result_ratio['dm_q30_bases_r2'])
    result_ratio['dm_soft_clipped_bases_diff'] = (result_ratio['dm_soft_clipped_bases_r1'] - result_ratio['dm_soft_clipped_bases_r2'])
    result_ratio['dm_reads_with_indel_diff'] = (result_ratio['dm_reads_with_indel_r1'] - result_ratio['dm_reads_with_indel_r2'])
    result_ratio['dm_mismatched_bases_diff'] = (result_ratio['dm_mismatched_bases_r1'] - result_ratio['dm_mismatched_bases_r2'])
    result_values['dm_mapped_bases'] = result_values['dm_mapped_bases_r1'] + result_values['dm_mapped_bases_r2']
    return (result_values, result_ratio)
        

def combine_rg_quality_stats(sample_readgroups, adapter_removals, adapters, merge_stats, dragmaps, dechimers):
    
    header = ['sample','readgroup']
    header_adapterr = ['ar_read_pairs', 'ar_aligned_fraction', 'ar_adapter_fraction', 'ar_retained_fraction', 'ar_average_read_length_retained']
    header_adapteri = ['ai_reads_aligned','ai_n_adapter','ai_adapter1','ai_adapter2']
    header_merge = ['merge_fragments', 'merge_alignments', 'merge_supplementary_alignments_ratio', 'merge_total_bp', 'merge_primary_soft_clipped_bp_ratio', 'merge_supplementary_bp_ratio',\
                    'merge_readded_fragments_ratio', 'merge_restored_read_ratio', 'merge_restored_bp_ratio']

    header_dm_value = ['dm_total_input_reads', 'dm_mapped_reads','dm_properly_paired_reads','dm_total_bases', 'dm_mapped_bases']
    header_dm_ratio = ['dm_mapped_reads', 'dm_unmapped_reads', 'dm_singleton_reads_itself_mapped_mate_unmapped',
        'dm_paired_reads_itself_mate_mapped', 'dm_not_properly_paired_reads_discordant', 'dm_paired_reads_mapped_to_different_chromosomes_mapq_10',
        'dm_reads_with_mapq_40_inf','dm_reads_with_mapq_30_40', 'dm_reads_with_mapq_20_30', 'dm_reads_with_mapq_10_20', 'dm_reads_with_mapq_0_10',
        'dm_reads_with_mapq_na_unmapped_reads', 'dm_reads_with_indel', 'dm_soft_clipped_bases', 'dm_mismatched_bases_excl_indels', 'dm_q30_bases', 'dm_soft_clipped_bases_diff', 
        'dm_reads_with_indel_diff', 'dm_mismatched_bases_diff', 'dm_q30_bases_diff']
    
    header_dechimer = ['dechimer_fragment_modified_ratio', 'dechimer_clip_ratio', 'dechimer_loose_end_clip_ratio', 'dechimer_all_sup_discarded_ratio', 'dechimer_min_alignment_length_unmap_ratio']

    
    header = header + header_adapterr + header_adapteri + header_merge + header_dm_value + header_dm_ratio + header_dechimer
    
    data = []

    for sample_rg, adapter_r, adapter_i, merge_stats, dragmap_stats, dechimer_stats in zip(sample_readgroups, adapter_removals, adapters, merge_stats, dragmaps, dechimers):
        sample, rg = sample_rg
        
        a = read_adapter_removal_stats(adapter_r)
        nrow = [sample,rg] + [a[h] for h in header_adapterr]
        
        a = read_adapter_identify_stats(adapter_i)
        nrow = nrow + [a[h] for h in header_adapteri]
        
        a = read_merge_stats(merge_stats)
        nrow = nrow + [a[h] for h in header_merge]
        
        a = read_dragmap_stats(dragmap_stats)
        nrow = nrow + [a[0][h] for h in header_dm_value] + [a[1][h] for h in header_dm_ratio]
                    
        a = read_dechimer_stats(dechimer_stats)
        nrow = nrow + [a[h] for h in header_dechimer]
                            
        data.append(nrow)
 
    return (tuple(header), data)


def combine_quality_stats(samples, genome_filenames, exome_filenames, vpca2, bam_extra_all, bam_extra_exome, pre_adapter, bait_bias, hsstats):
    header = ['sample', 'total_sequences', 'avg_quality', 'average_length', 'max_length', 'error_rate',\
              'insert_size_avg', 'insert_size_std',\
              'reads_unmapped', 'reads_properly_paired', 'reads_mq0', 'reads_duplicated',\
              'pairs_inward', 'pairs_outward', 'pairs_other', 'pairs_diff_chromosomes',\
              'bases_total', 'bases_mapped',\
              'qual_by_cycle','base_fractions', 'insert_maxsize', 'pair_inward_maxsize_frac','pair_outward_maxsize_frac','pair_other_maxsize_frac',\
              'pair_outward_size0_frac','pair_other_size0_frac','gc35_depth','gc50_depth',\
              'exome_total_sequences', 'exome_reads_mq0', 'exome_insert_size_avg', 'exome_insert_size_std','exome_pairs_inward', 'exome_pairs_outward', 'exome_pairs_other', 'exome_pairs_diff_chromosomes',\
              'exome_insertions_frac','exome_deletions_frac','exome_insertions_length_avg','exome_deletions_length_avg','exome_coverage1_frac','exome_coverage1000plus_frac','exome_gc35_depth','exome_gc50_depth',\
              'pca2_freemix',\
              'unmapped_ratio_all','mqual20_ratio_all','secondary_ratio_all','supplementary_ratio_all','sup_diffchrom_ratio_all','duplicate_ratio_all','duplicate_supplement_ratio_all','soft_clipped_bp_ratio_all',\
              'aligned_bp_ratio_all','inserted_bp_ratio_all','deleted_bp_ratio_all','total_bp_all','soft_clipped_bp_ratio_filter_50_all','soft_clipped_bp_ratio_filter_60_all','soft_clipped_bp_ratio_filter_70_all',\
              'aligned_bp_ratio_filter_50_all','aligned_bp_ratio_filter_60_all','aligned_bp_ratio_filter_70_all','poly_a_all','poly_g_all', 'illumina_adapter_all','pcr_adapter_1_all','pcr_adapter_2_all','nextera_all',\
              'unmapped_ratio_exome','mqual20_ratio_exome','secondary_ratio_exome','supplementary_ratio_exome','sup_diffchrom_ratio_exome','duplicate_ratio_exome','duplicate_supplement_ratio_exome',\
              'soft_clipped_bp_ratio_exome','aligned_bp_ratio_exome','inserted_bp_ratio_exome','deleted_bp_ratio_exome','total_bp_exome','soft_clipped_bp_ratio_filter_50_exome','soft_clipped_bp_ratio_filter_60_exome',\
              'soft_clipped_bp_ratio_filter_70_exome','aligned_bp_ratio_filter_50_exome','aligned_bp_ratio_filter_60_exome','aligned_bp_ratio_filter_70_exome','poly_a_exome','poly_g_exome', 'illumina_adapter_exome','pcr_adapter_1_exome',\
              'pcr_adapter_2_exome','nextera_exome', \
              'pre_adapter_ac', 'pre_adapter_ag','pre_adapter_at', 'pre_adapter_ca', 'pre_adapter_cg', 'pre_adapter_ct', 'pre_adapter_ga', 'pre_adapter_gc', 'pre_adapter_gt', 'pre_adapter_ta', 'pre_adapter_tc', 'pre_adapter_tg',\
              'bait_bias_ac', 'bait_bias_ag','bait_bias_at', 'bait_bias_ca', 'bait_bias_cg', 'bait_bias_ct', 'bait_bias_ga', 'bait_bias_gc', 'bait_bias_gt', 'bait_bias_ta', 'bait_bias_tc', 'bait_bias_tg']
              
    bamstats_header = ['unmapped_ratio','mqual20_ratio','secondary_ratio','supplementary_ratio','sup_diffchrom_ratio','duplicate_ratio','duplicate_supplement_ratio','soft_clipped_bp_ratio',\
                        'aligned_bp_ratio','inserted_bp_ratio','deleted_bp_ratio','total_bp','soft_clipped_bp_ratio_filter_50','soft_clipped_bp_ratio_filter_60','soft_clipped_bp_ratio_filter_70',\
                        'aligned_bp_ratio_filter_50','aligned_bp_ratio_filter_60','aligned_bp_ratio_filter_70','poly_a','poly_g', 'illumina_adapter','pcr_adapter_1','pcr_adapter_2','nextera']    
    
    hsstat_header = ['bait_design_efficiency',
                    'on_bait_bases','near_bait_bases','off_bait_bases',
                    'pct_off_bait',	'mean_bait_coverage','pct_usable_bases_on_bait','pct_usable_bases_on_target',
                    'fold_enrichment',
                    'on_target_bases','mean_target_coverage',	'median_target_coverage','zero_cvg_targets_pct',
                    'pct_exc_dupe',	'pct_exc_mapq',	'pct_exc_baseq','pct_exc_overlap','pct_exc_off_target',
                    'pct_target_bases_1x',	'pct_target_bases_2x',	'pct_target_bases_10x',	'pct_target_bases_20x',
                    'pct_target_bases_30x',	'pct_target_bases_50x',	'pct_target_bases_100x','pct_target_bases_250x',
                    'pct_target_bases_500x',
                    'at_dropout','gc_dropout',
                    'het_snp_sensitivity']
    header = header +  hsstat_header

    data = []

    for sample, genome_filename, exome_filename, v2,ba,be, pre_ad, bait_b, hsstat in zip(samples, genome_filenames, exome_filenames, vpca2, bam_extra_all, bam_extra_exome, pre_adapter, bait_bias, hsstats):
        print('Reading statistics for sample: ' + sample)
        stats = read_stats_full(genome_filename)
        
        nrow = [sample, stats['raw total sequences'], stats['average quality'], stats['average length'], stats['maximum length'], stats['error rate'],\
                        stats['insert size average'], stats['insert size standard deviation'],\
                        stats['reads unmapped'], stats['reads properly paired'], stats['reads MQ0'], stats['reads duplicated'],\
                        stats['inward oriented pairs'], stats['outward oriented pairs'], stats['pairs with other orientation'], stats['pairs on different chromosomes'],\
                        stats['total length'], stats['bases mapped (cigar)'], \
                        stats['qual_by_cycle'], stats['base_fractions'], stats['insert_maxsize'],\
                        stats['pair_inward_maxsize'], stats['pair_outward_maxsize'], stats['pair_other_maxsize'], stats['pair_outward_size0'], stats['pair_other_size0'],\
                        stats['gc_depth35'], stats['gc_depth50']]

        exome_stats = read_stats_full(exome_filename)
        nrow = nrow + [exome_stats['raw total sequences'], exome_stats['reads MQ0'], exome_stats['insert size average'], exome_stats['insert size standard deviation'],\
                       exome_stats['inward oriented pairs'], exome_stats['outward oriented pairs'], exome_stats['pairs with other orientation'], exome_stats['pairs on different chromosomes'],\
                       exome_stats['insertions_frac'], exome_stats['deletions_frac'], exome_stats['insertions_avg'], exome_stats['deletions_avg'], exome_stats['coverage1'], exome_stats['coverage1000'],\
                       exome_stats['gc_depth35'], exome_stats['gc_depth50']]
        
        #nrow = nrow + [read_contam(veurbam_filename)]
        #nrow = nrow + [read_contam(v1000gbam_filename)]
        # where v2 is going from?
        nrow = nrow + [read_contam(v2)]
        #nrow = nrow + [read_contam(v4)]
        sdata = read_bamstat(ba)
        nrow = nrow + [sdata[h] for h in bamstats_header]
        sdata = read_bamstat(be)
        nrow = nrow + [sdata[h] for h in bamstats_header]
        
        nrow = nrow + read_artifacts(pre_ad)
        nrow = nrow + read_artifacts(bait_b)

        hsstat = read_hsstats(hsstat)[0]
        nrow = nrow + [hsstat[h] for h in hsstat_header]

        data.append(tuple(nrow))
    return (tuple(header), data)

def read_bamstat(filename):
    print('BAMSTAT', filename)
    
    with open(filename, 'r') as f:
        header = f.readline()
        data = f.readline()
        fields = header.strip().lstrip('#').split('\t')
        data = data.strip().split('\t')
        res = dict(zip(fields, data))
    print(res)        
    return res


    
def correct_between(*datas):
    pass
    return (datas, None)




def combine_oxo_stats(samples, pre_adapter, bait_bias):

    header = ['sample']

    contexts = ['aa','ac','ag','at','ca','cc','cg','ct','ga','gc','gg','gt', 'ta','tc','tg','tt']

    for i in ['A','C','G','T']:
        for j in ['A','C','G','T']:
            if i != j:
                for con in contexts:
                    header.append('pre_%s_%s_%s' % (i, j, con))

    for i in ['A','C','G','T']:
        for j in ['A','C','G','T']:
            if i != j:
                for con in contexts:
                    header.append('bait_%s_%s_%s' % (i, j, con))
                
    data = []

    for sample, pre_ad, bait_b in zip(samples, pre_adapter, bait_bias):
        print('Reading statistics for sample: ' + sample)
        nrow = [sample]

        preadapter = read_detail_artifacts(pre_ad)
        bait = read_detail_artifacts(bait_b)


        for i in ['a','c','g','t']:
            for j in ['a','c','g','t']:
                if i != j:
                    for con in contexts:
                        fr,fa, rr, ra  = preadapter[(i,j,con)]
                        e1 = float(fa)/float(fa + fr) 
                        e2 = float(ra)/float(ra + rr) 
                        pe1 = numpy.log10(e1) * -10.0
                        pe2 = numpy.log10(e2) * -10.0

                        if fa > ra:
                            e3 = (float(fa) - float(ra)) / (float(fa) + float(ra) + float(fr) + float(rr))
                            pe3 = min(numpy.log10(e3) * -10.0, 100.0)
                        else:
                            pe3 = 100.0
                        nrow.append('%.3f;%.3f;%.3f' % (pe1,pe2, pe3))

        for i in ['a','c','g','t']:
            for j in ['a','c','g','t']:
                if i != j:
                    for con in contexts:
                        fr,fa,rr,ra = bait[(i,j,con)]
                        e1 = float(fa)/float(fa + fr) 
                        e2 = float(ra)/float(ra + rr) 
                        pe1 = numpy.log10(e1) * -10.0
                        pe2 = numpy.log10(e2) * -10.0
                        if e1 > e2:
                            pe3 = min(numpy.log10((e1 - e2)) * -10.0,100.0)
                        else:
                            pe3 = 100.0
                        nrow.append('%.3f;%.3f;%.3f' % (pe1,pe2,pe3))
        data.append(tuple(nrow))
    return (tuple(header), data)

