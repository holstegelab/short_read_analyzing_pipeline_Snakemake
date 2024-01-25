import csv
import sys


#load bed file
bed_file = sys.argv[1]
bed = open(bed_file, 'r')
bed_reader = csv.reader(bed, delimiter='\t')


#load genome file (chromosome lengths)
genome_file = sys.argv[2]
genome = open(genome_file, 'r')
genome_reader = csv.reader(genome, delimiter='\t')
genome_dict = { row[0]:int(row[1]) for row in genome_reader }

def slop(bed_record, genome_dict, add_start=0, add_stop=0):
    #now add 1000 bp to each side of the bed record
    #if the bed record is within 1000 bp of the start or end of the chromosome, then just extend to the end of the chromosome
    
    #get the chromosome
    chrom = bed_record[0]
    #get the start and end positions
    start = int(bed_record[1])
    end = int(bed_record[2])
    
    #get the chromosome length
    chrom_length = genome_dict[chrom]


    start = max(0, start-add_start)
    end = min(chrom_length-1, end+add_stop)
    return (chrom, start, end)


bedrecords = [(row[0], int(row[1]), int(row[2])) for row in bed_reader]
bedrecords[0] = slop(bedrecords[0], genome_dict, add_start=1000, add_stop=0)
bedrecords[-1] = slop(bedrecords[-1], genome_dict, add_start=0, add_stop=1000)

#write the bed records to a new file
new_filename = '.'.join(bed_file.split('.')[:-1]) + '.padded.bed'
with open(new_filename, 'w') as outfile:
    for record in bedrecords:
        outfile.write('\t'.join([str(x) for x in record]) + '\n')



