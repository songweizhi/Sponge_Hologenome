from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def get_Porifera_seqs(seq_in, seq_out):
    taxon_higher_Porifera = 'Eukaryota;Amorphea;Obazoa;Opisthokonta;Holozoa;Choanozoa;Metazoa;Porifera;'
    seq_out_handle = open(seq_out, 'w')
    for each_seq in SeqIO.parse(seq_in, 'fasta'):
        seq_desc = each_seq.description
        if taxon_higher_Porifera in seq_desc:
            seq_out_handle.write('>%s\n' % seq_desc.replace(taxon_higher_Porifera, ''))
            seq_out_handle.write('%s\n' % str(each_seq.seq).replace('U', 'T'))
    seq_out_handle.close()


def plot_seq_len_distribution(seq_in, plot_out):
    length_list = []
    for each_seq in SeqIO.parse(seq_in, 'fasta'):
        length_list.append(len(each_seq.seq))
    plt.hist(length_list, bins=100)
    plt.savefig(plot_out, dpi=300)
    plt.close()


def export_seq_description(seq_in, seq_des_txt):
    seq_des_txt_handle = open(seq_des_txt, 'w')
    for each_seq in SeqIO.parse(seq_in, 'fasta'):
        seq_id = each_seq.id
        seq_desc = each_seq.description[len(seq_id)+1:]
        seq_des_txt_handle.write('%s\t%s\n' % (seq_id, seq_desc))
    seq_des_txt_handle.close()


###################################################### file in/out #####################################################

silva_wd = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_SILVA'

# file in SSU
SSUParc_seq                   = '%s/SILVA_138.1_SSUParc_tax_silva.fasta'                                      % silva_wd
SSURef_NR99_seq               = '%s/SILVA_138.1_SSURef_NR99_tax_silva.fasta'                                  % silva_wd

# file in LSU
LSUParc_seq                   = '%s/SILVA_138.1_LSUParc_tax_silva.fasta'                                      % silva_wd
LSURef_NR99_seq               = '%s/SILVA_138.1_LSURef_NR99_tax_silva.fasta'                                  % silva_wd

# file out SSU
SSUParc_Sponge_18S_seq        = '%s/SILVA_Porifera/SILVA_138.1_SSUParc_tax_silva_Sponge_18S.fasta'            % silva_wd
SSUParc_Sponge_18S_tax        = '%s/SILVA_Porifera/SILVA_138.1_SSUParc_tax_silva_Sponge_18S_taxon.txt'        % silva_wd
SSUParc_Sponge_18S_png        = '%s/SILVA_Porifera/SILVA_138.1_SSUParc_tax_silva_Sponge_18S.png'              % silva_wd
SSURef_NR99_Sponge_18S_seq    = '%s/SILVA_Porifera/SILVA_138.1_SSURef_NR99_tax_silva_Sponge_18S.fasta'        % silva_wd
SSURef_NR99_Sponge_18S_tax    = '%s/SILVA_Porifera/SILVA_138.1_SSURef_NR99_tax_silva_Sponge_18S_taxon.txt'    % silva_wd
SSURef_NR99_Sponge_18S_png    = '%s/SILVA_Porifera/SILVA_138.1_SSURef_NR99_tax_silva_Sponge_18S.png'          % silva_wd

# file out LSU
LSUParc_Sponge_28S_seq        = '%s/SILVA_Porifera/SILVA_138.1_LSUParc_tax_silva_Sponge_28S.fasta'            % silva_wd
LSUParc_Sponge_28S_tax        = '%s/SILVA_Porifera/SILVA_138.1_LSUParc_tax_silva_Sponge_28S_taxon.txt'        % silva_wd
LSUParc_Sponge_28S_png        = '%s/SILVA_Porifera/SILVA_138.1_LSUParc_tax_silva_Sponge_28S.png'              % silva_wd

LSURef_NR99_Sponge_28S_seq    = '%s/SILVA_Porifera/SILVA_138.1_LSURef_NR99_tax_silva_Sponge_28S.fasta'        % silva_wd
LSURef_NR99_Sponge_28S_tax    = '%s/SILVA_Porifera/SILVA_138.1_LSURef_NR99_tax_silva_Sponge_28S_taxon.txt'    % silva_wd
LSURef_NR99_Sponge_28S_png    = '%s/SILVA_Porifera/SILVA_138.1_LSURef_NR99_tax_silva_Sponge_28S.png'          % silva_wd


################################################ get Porifera sequences ################################################

# get_Porifera_seqs(SSUParc_seq, SSUParc_Sponge_18S_seq)
# get_Porifera_seqs(SSURef_NR99_seq, SSURef_NR99_Sponge_18S_seq)

# get_Porifera_seqs(LSUParc_seq, LSUParc_Sponge_28S_seq)
# get_Porifera_seqs(LSURef_NR99_seq, LSURef_NR99_Sponge_28S_seq)

##################################### get length distribution of Porifera sequences ####################################

# plot_seq_len_distribution(SSUParc_Sponge_18S_seq, SSUParc_Sponge_18S_png)
# plot_seq_len_distribution(SSURef_NR99_Sponge_18S_seq, SSURef_NR99_Sponge_18S_png)

# plot_seq_len_distribution(LSUParc_Sponge_28S_seq, LSUParc_Sponge_28S_png)
# plot_seq_len_distribution(LSURef_NR99_Sponge_28S_seq, LSURef_NR99_Sponge_28S_png)

############################################## write out seq description ###############################################

# export_seq_description(SSUParc_Sponge_18S_seq, SSUParc_Sponge_18S_tax)
# export_seq_description(LSUParc_Sponge_28S_seq, LSUParc_Sponge_28S_tax)

# export_seq_description(SSURef_NR99_Sponge_18S_seq, SSURef_NR99_Sponge_18S_tax)
# export_seq_description(LSURef_NR99_Sponge_28S_seq, LSURef_NR99_Sponge_28S_tax)

########################################################################################################################

