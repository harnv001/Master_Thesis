#!/usr/bin/env python

"""
Program that extract the domain information from plantiSMASH output (genecluster.js) and map to the sequence list
to construct the protein domain format (based on http://itol.embl.de/help.cgi#domains)
"""

from __future__ import division

__author__ = "Yosapol Harnvanichvech (950416798110)"
__date__   = "09/06/2016"
__email__  = "yosapol.harnvanichvech@wur.nl"

#Imports
import json
import argparse
import subprocess
from database import database

#Functions
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sequence_name',
                        dest = 'sequence_name',
                        required = True,
                        type = str,
                        nargs = 1,
                        help = "put sequence_name in the textfile with one name per line ex. LOCUS|NAME-plant name")
    parser.add_argument('-json_directory',
                        dest ='json_directory',
                        required = True,
                        nargs = '+',
                        help = "put directoy of json file for extract domain information, based on sequence name")
    parser.add_argument('-output_name',
                        dest='output_name',
                        required = True,
                        type = str,
                        nargs = 1,
                        help = "put the file name of output with directory")
    parser.add_argument('-heatmap_name',
                        dest = 'heatmap_name',
                        required = True,
                        type = str,
                        nargs=1,
                        help ="put the heatmap_name")
    return parser.parse_args()

def extract_sequence(sequence_name):
    """extract the name of sequence and keep in the list.

          Keyword arguments:
              sequence_name -- sequence_name in the textfile with one name per line"

          Returns:
              sequence_list -- sequence_name in the list
      """
    sequence_list = []
    for name in sequence_name:
        if name.startswith(">"):
            name = name.replace(">","").strip("/n")
            sequence_list +=[name]
    # print sequence_list
    return sequence_list

def json_trimheader(genecluster_json):

    """Trim header and end line of json file from plantiSMASH output (genecluster.js).

        Keyword arguments:
            genecluster_json -- directory of genecluster.js

        Returns:
            json_output -- directory of trimmed header and end-line in json format
    """

    json_folder = []

    for folder in genecluster_json:
        path_output = "/".join(folder.split("/")[:-1]) +"/"
        trim_file = open(path_output+ "genecluster_trim.js", "w")
        json_folder += [path_output+"genecluster_trim.js"]
        with open(folder, "r") as genecluster_file:
            for line in genecluster_file:
                line = line.strip()
                line = line.strip("var geneclusters =")
                line = line.strip("var details_data = }")
                if line.endswith(";"):
                    line = line.replace(";","").strip()
                if line.endswith("{}"):
                    line = line.replace("{}", "").strip()
                trim_file.write(str(line))

    return json_folder

def extract_location(json_folder):
    """ check the location  of each sequence and keep in dictionary

        Keyword arguments:
            json_folder -- directory of trimmed header and end-line in json format (genecluster_trim.js.)

        Returns:
            domain_data -- dictionary of domain data using name as a keys and domain information as a value.
    """
    domain_data = {}
    for folder in json_folder:
        json_trimfile = json.loads(open(folder).read())
        all_cluster = json_trimfile.keys() # number of cluster in sequences
        # database = {"p450": ["Cytochrome_450", "#ff0000", "RE"],
        #         "Terpene_synth_C": ["Terpene synthase", "#35F800", "RE"],
        #         "SQHop_cyclase_C": ["Terpene synthase", "#35F800", "RE"],
        #         "SQHop_cyclase_N": ["Terpene synthase", "#35F800", "RE"],
        #         "SQHop_cyclase_C,SQHop_cyclase_N": ["Terpene synthase", "#35F800", "RE"],
        #         "Terpene_synth": ["Terpene synthase", "#35F800", "RE"],
        #         "Lycopene_cycl": ["Terpene synthase", "#35F800", "RE"],
        #         "Cu_amine_oxid": ["Copper amine oxidase", "#E9C63E", "RE"],
        #         "Bet_v_1": ["Pictet-Spengler enzyme (Bet v1)", "#735100", "RE"],
        #         "Glycos_transf_1": ["Glycosyltransferase", "#F477A6", "RE"],
        #         "Glycos_transf_2": ["Glycosyltransferase", "#F477A6", "RE"],
        #         "Glyco_transf_28": ["Glycosyltransferase", "#F477A6", "RE"],
        #         "UDPGT": ["Glycosyltransferase", "#F477A6", "RE"],
        #         "UDPGT_2": ["Glycosyltransferase", "#F477A6", "RE"],
        #         "Chal_sti_synt_C": ["Ketosynthase", "#75FFD8", "RE"],
        #         "Chal_sti_synt_N": ["Ketosynthase", "#75FFD8", "RE"],
        #         "SE": ["Squalene epoxidase", "#206B14", "RE"], "COesterase": ["COesterase", "#A221FF", "RE"],
        #         "Methyltransf_2": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "Methyltransf_3": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "Methyltransf_7": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "Methyltransf_11": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "cMT": ["Methyltransferase", "#CF0FFF", "RE"], "nMT": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "oMT": ["Methyltransferase", "#CF0FFF", "RE"],
        #         "Transferase": ["BAHD acyltransferase", "#0003EC", "RE"],
        #         "Peptidase_S10": ["Scl acyltransferase", "#0c126d", "RE"],
        #         "Epimerase": ["Epimerase", "#F245BD", "RE"],
        #         "adh_short": ["Oxidoreductase", "#4DD0E1", "RE"],
        #         "adh_short_C2": ["Oxidoreductase", "#4DD0E1", "RE"],
        #         "NAD_binding_1": ["Oxidoreductase", "#4DD0E1", "RE"],
        #         "GMC_oxred_N": ["Oxidoreductase", "#4DD0E1", "RE"],
        #         "GMC_oxred_C": ["Oxidoreductase", "#4DD0E1 ", "RE"], "DIOX_N": ["Dioxygenase", "#F5B2FF", "RE"],
        #         "2OG-FeII_Oxy": ["Dioxygenase", "#F5B2FF", "RE"],
        #         "AMP-binding": ["CoA-ligase", "#8F400B", "RE"], "Amino_oxidase": ["Amino oxidase", "#2C113A", "RE"],
        #         "Aminotran_1_2": ["Aminotransferase", "#E42E00", "RE"],
        #         "Aminotran_3": ["Aminotransferase", "#E42E00", "RE"],
        #         "Prenyltrans": ["Prenyltransferase", "#552288", "RE"],
        #         "Prenyltransf": ["Prenyltransferase", "#552288", "RE"],
        #         "UbiA": ["Prenyltransferase", "#552288", "RE"],
        #         "Str_synth": ["Strictosidine synthase-like", "#dd71dd", "RE"],
        #         "PRISE": ["PRISE enzymes", "#13C600", "RE"], "Dirigent": ["Dirigent enzymes", "#606000", "RE"]
        #         }
        # extract all locus tag with location(start,stop)
        for cluster in all_cluster:
            locus_tmp = []
            data_tmp = []
            name_cluster = str(cluster)
            location_start = json_trimfile[cluster]['start']
            location_end = json_trimfile[cluster]['end']
            seq_length = int(location_end - location_start)
            for gene in json_trimfile[cluster]['orfs']:
                if "Signature pHMM" in gene['description']:
                    locus_tag = gene['locus_tag']
                    gene_start = gene['start'] # position gene start
                    gene_end = gene['end'] # position gene end
                    scale_start = gene_start - location_start # adjust scale of figure to the same length of domain
                    scale_end = gene_end - location_start  # adjust scale of figure to the same length of domain
                    #end = str(gene['end'])
                    domain = gene['domains']
                    domain = [domain_name.replace("plants/", "") for domain_name in domain] # trim "plant/"
                    for domain_name in [domain[0]]:
                        try:
                            domain_figure = database[domain_name]
                        except KeyError:
                            domain_figure = [str(domain_name), "#800000", "RE"] #  others biosynthetic genes
                        name_domain = ",".join(str(name) for name in domain)
                    locus_tmp.append(str(locus_tag))

                    data_tmp.append([seq_length,
                                     gene_start,
                                     gene_end,
                                     scale_start,
                                     scale_end,
                                     name_domain,
                                     name_cluster] + domain_figure)
            for locus in locus_tmp:
                domain_data[locus] = data_tmp
    for folder in json_folder:
        cmd = "rm {}".format(folder)
        if subprocess.check_call(cmd, shell=True) == 0:
            print ("delete"+ str(folder))

    return domain_data

def domain_format(sequence_name,output_name,domain_data,heatmap_name):
    """ If sequence name were contained in domain datam write protein domain format based on protein domain 
        format in (http://itol.embl.de/help.cgi#domains)
        
        Keyword arguments:
            domain_data -- dictionary of domain data using name as a keys and domain information as a value
            domain_data -- dictionary of domain data using name as a keys and domain information as a value.
            output_name -- directory of outputfile.txt
                
        Returns:
            file output.txt in output_name directory
    """
# write header format of protein_format file
    protein_format = open(output_name[0], "w")
    protein_format.write("DATASET_DOMAINS" + '\n')
    protein_format.write("SEPARATOR TAB" + '\n')
    protein_format.write("BACKBONE_COLOR" + '\t' + "#000000" + '\n')
    protein_format.write("BACKBONE_HEIGHT" + '\t' + "3" + '\n')
    protein_format.write("BORDER_WIDTH" + '\t' + "0" + '\n')
    protein_format.write("COLOR" + '\t' + "#bebada" + '\n')
    protein_format.write("DATASET_LABEL" + '\t' + heatmap_name[0] + '\n')
    protein_format.write("HEIGHT_FACTOR" + '\t' + "1" + '\n')
    protein_format.write("DATA" + '\n')

    name_domain = open(sequence_name[0], "r")
    for name in name_domain:
        locus_name = name
        key_name = name.partition('|')[0] # The format of gene tag in this script is "id|name of genome file-plant name"
        write_tmp = ""
        if key_name in domain_data.keys():
            protein_inf = domain_data[key_name]
            for protein_data in protein_inf:
                seq_length = str(protein_data[0])
                seq_start = str(protein_data[3])
                seq_end = str(protein_data[4])
                label = str(protein_data[7])
                color = str(protein_data[8])
                shape = str(protein_data[9])
                write_tmp += shape + "|" + seq_start + "|" + seq_end + "|" + color + "|" + label + '\t'
            protein_format.write("\t".join([locus_name.strip("\n"), seq_length, write_tmp]) + "\n")
    print "finish writing protein format"


if __name__ == "__main__":

    #Defining the fileName variables
    global arguments
    arguments = get_arguments()
    sequence_name = arguments.sequence_name
    genecluster_json = arguments.json_directory
    output_name = arguments.output_name
    heatmap_name = arguments.heatmap_name
    json_folder = json_trimheader(genecluster_json)
    domain_data = extract_location(json_folder)
    domain_format(sequence_name, output_name, domain_data, heatmap_name)

