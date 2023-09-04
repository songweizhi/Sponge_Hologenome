import os
import glob
from Bio import SeqIO
import numpy as np

def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def calculate_coding_density(gff_file):

    '''
    Calculates the genome density for the given GFF
    '''

    # read GFF
    lines  = [line.strip().split("\t") for line in open(gff_file).readlines()]
    #print(lines[1])

    # get length of the genome
    length = float(lines[1][0].split(";seqlen=")[1].split(";")[0])

    # store start and ends of the ORFs
    starts = [int(line[3]) for line in lines if not line[0].startswith("#")]
    ends   = [int(line[4]) for line in lines if not line[0].startswith("#")]

    # compute average length
    lengths = list()
    for start, end in zip(starts, ends):
        lengths.append(int(end)-int(start))
    mean_len = "{:.3f}".format(np.mean(lengths))

    # compute coding length. It collapses overlaping ORFs in the same
    # interval
    intervals = [[s,e] for s,e in zip(starts, ends)]

    intervals.sort(key=lambda interval: interval[0])
    merged = [intervals[0]]
    for current in intervals:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)

    # coding length
    coding_bases = float(0)
    for interval in merged:
        coding_bases += 1 + (interval[1] - interval[0])

    # coding density
    density = "{:.3f}".format(coding_bases/length)

    # return id, 11/TAG/TGA, seqlen, density, n_prots, prots_len
    return [str(int(length)), density, str(len(starts)), mean_len]


########################################################################################################################

# file in
gnm_dir = '/Users/songweizhi/Desktop/Sponge_2023_08_25/Nitrosopumilaceae_50_5_dRep97_195'
gnm_ext = 'fna'
ffn_dir = '/Users/songweizhi/Desktop/Sponge_2023_08_25/Nitrosopumilaceae_50_5_dRep97_195_ffn'

# file out
coding_density_txt = '/Users/songweizhi/Desktop/Sponge_2023_08_25/Nitrosopumilaceae_50_5_dRep97_195_coding_density.txt'


'''

cd /Users/songweizhi/Desktop/Sponge_2023_08_25
BioSAK get_gnm_size -i Nitrosopumilaceae_50_5_dRep97_195 -x fna > Nitrosopumilaceae_50_5_dRep97_195_gnm_size.txt
BioSAK iTOL -SimpleBar -lv Nitrosopumilaceae_50_5_dRep97_195_gnm_size.txt -scale 0-1-2-3 -lt Size -o Nitrosopumilaceae_50_5_dRep97_195_gnm_size_iTOL.txt
BioSAK iTOL -SimpleBar -lv Nitrosopumilaceae_50_5_dRep97_195_gnm_size_minus_0.5.txt -scale 0.5-1-1.5-2 -lt Size -o Nitrosopumilaceae_50_5_dRep97_195_gnm_size_minus_0.5_iTOL.txt

cd /Users/songweizhi/Desktop/Sponge_2023_08_25
BioSAK iTOL -SimpleBar -lv Nitrosopumilaceae_50_5_dRep97_195_coding_density.txt -scale 0-10-20-30 -lt Coding_Density -o Nitrosopumilaceae_50_5_dRep97_195_coding_density_iTOL.txt
BioSAK iTOL -SimpleBar -lv Nitrosopumilaceae_50_5_dRep97_195_coding_density_minus_70.txt -scale 0-10-20-30 -lt Coding_Density -o Nitrosopumilaceae_50_5_dRep97_195_coding_density_minus_70_iTOL.txt

'''

########################################################################################################################

gnm_file_re = '%s/*.%s' % (gnm_dir, gnm_ext)
gnm_file_list = glob.glob(gnm_file_re)

coding_density_txt_handle = open(coding_density_txt, 'w')
for each_gnm in gnm_file_list:

    _, f_base, _ = sep_path_basename_ext(each_gnm)
    ffn_file = '%s/%s.ffn' % (ffn_dir, f_base)

    gnm_size_bp = 0
    for each_seq in SeqIO.parse(each_gnm, 'fasta'):
        gnm_size_bp += len(each_seq.seq)

    orf_len_bp = 0
    for each_gene in SeqIO.parse(ffn_file, 'fasta'):
        orf_len_bp += len(each_gene.seq)

    coding_seq_pct = orf_len_bp*100/gnm_size_bp
    coding_seq_pct = float("{0:.2f}".format(coding_seq_pct))
    coding_density_txt_handle.write('%s\t%s\n' % (f_base, coding_seq_pct))

