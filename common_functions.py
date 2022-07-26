
def gtdb_tax_str_parser(gtdb_tax_str, taxon_rank):
    gtdb_tax_str_split = gtdb_tax_str.split(';')
    taxon_to_return = ''
    if taxon_rank == 'd':
        taxon_to_return = gtdb_tax_str_split[0]
    elif taxon_rank == 'p':
        taxon_to_return = gtdb_tax_str_split[1]
    elif taxon_rank == 'c':
        taxon_to_return = gtdb_tax_str_split[2]
    elif taxon_rank == 'o':
        taxon_to_return = gtdb_tax_str_split[3]
    elif taxon_rank == 'f':
        taxon_to_return = gtdb_tax_str_split[4]
    elif taxon_rank == 'g':
        taxon_to_return = gtdb_tax_str_split[5]
    elif taxon_rank == 's':
        taxon_to_return = gtdb_tax_str_split[6]
    return taxon_to_return
