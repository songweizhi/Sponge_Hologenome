from Bio import SeqIO


msa_file        = '/Users/songweizhi/Desktop/Sponge_2023_08_25/top25.concatenated.phy.fasta'
msa_seq_pct_txt = '/Users/songweizhi/Desktop/Sponge_2023_08_25/top25.concatenated.phy.fasta.gap.txt'


msa_seq_pct_txt_handle = open(msa_seq_pct_txt, 'w')
for each_seq in SeqIO.parse(msa_file,'fasta'):
    msa_seq = each_seq.seq
    dash_num = msa_seq.count('-')
    seq_pct = 100 - dash_num*100/len(msa_seq)
    seq_pct = float("{0:.2f}".format(seq_pct))
    print('%s\t%s' % (each_seq.id, seq_pct))
    msa_seq_pct_txt_handle.write('%s\t%s\n' % (each_seq.id, seq_pct))
msa_seq_pct_txt_handle.close()

'''

cd /Users/songweizhi/Desktop/Sponge_2023_08_25
BioSAK iTOL -SimpleBar -lv top25.concatenated.phy.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top25.concatenated.phy.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top25.concatenated.phy.fasta.RmHeteSite_25.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top25.concatenated.phy.fasta.RmHeteSite_25.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top25.concatenated.phy.fasta.RmHeteSite_50.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top25.concatenated.phy.fasta.RmHeteSite_50.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top25.concatenated.phy.fasta.RmHeteSite_75.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top25.concatenated.phy.fasta.RmHeteSite_75.fasta.gap.iTOL.txt

BioSAK iTOL -SimpleBar -lv top50.concatenated.phy.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top50.concatenated.phy.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top50.concatenated.phy.fasta.RmHeteSite_25.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top50.concatenated.phy.fasta.RmHeteSite_25.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top50.concatenated.phy.fasta.RmHeteSite_50.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top50.concatenated.phy.fasta.RmHeteSite_50.fasta.gap.iTOL.txt
BioSAK iTOL -SimpleBar -lv top50.concatenated.phy.fasta.RmHeteSite_75.fasta.gap.txt -scale 0-25-50-75-100 -lt Gap -o top50.concatenated.phy.fasta.RmHeteSite_75.fasta.gap.iTOL.txt

'''