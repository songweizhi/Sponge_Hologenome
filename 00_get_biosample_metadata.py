import os
import argparse
import multiprocessing as mp


def get_biosample_metadata(args):

    file_in = args['i']

    cmd_list = []
    for each_id in open(file_in):

        biosample_id = each_id.strip()
        cmd = 'esearch -db biosample -query %s | esummary | xtract -pattern DocumentSummary -element Accession -element Date -first Title -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute -group Attribute -if Attribute@harmonized_name -equals "host" -element Attribute > %s_metadata.txt' % (biosample_id, biosample_id)
        cmd = 'esearch -db biosample -query %s -api_key fcbcd87b5c19888293ae07ca6c8da36d5108 | esummary -api_key fcbcd87b5c19888293ae07ca6c8da36d5108 | xtract -pattern DocumentSummary -element Accession -first Title -group Attribute -if Attribute@harmonized_name -equals "host" -element Attribute -group Attribute -if Attribute@harmonized_name -equals "isolation_source" -element Attribute > %s_metadata.txt' % (biosample_id, biosample_id)
        cmd_list.append(cmd)

    pool = mp.Pool(processes=12)
    pool.map(os.system, cmd_list)
    pool.close()
    pool.join()


if __name__ == '__main__':

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', required=True, help='biosample id file, one id per line')
    args = vars(arg_parser.parse_args())
    get_biosample_metadata(args)

'''
split -l 2500 combined_ar53_bac120_metadata_r207_col_biosample_uniq.tsv subset
'''