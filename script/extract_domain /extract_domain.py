#!/usr/bin/env python

"""
Program that extract the specific domain genecluster.js (plantiSMASH result)

"""
from __future__ import division

__author__ = "Yosapol Harnvanichvech (950416798110)"
__date__   = "09/06/2016"
__email__  = "yosapol.harnvanichvech@wur.nl"
"align_domain.py were modified from Satria Kautsar (PhD., Bioinformatics group, WUR)"

#Imports
import json
import re
import argparse
from align_domain import main

#Functions
def get_arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument('-input',
                        dest = 'input',
                        type = str,
                        required = True,
                        nargs = '+',
                        help = "put multiple directory of genecluster.js, using space as a seperator")
    parser.add_argument('-fasta_output',
                        dest = 'fasta_output',
                        type = str,
                        required = True,
                        help = "put directory of output")
    parser.add_argument('-domain',
                        dest = 'domain',
                        type = str,
                        required = True,
                        help = "put the detected domain")
    parser.add_argument('-known_domain',
                        dest = 'known_domain',
                        nargs='+',
                        type = str,
                        default=[],
                        help = "put known-domain.hmm, and seperated by space")

    return parser.parse_args()

def json_trimheader(genecluster_json):

    """Trim header and end line of json file from plantiSMASH output (genecluster.js).

        Keyword arguments:
            genecluster_json -- directory of genecluster.js

        Returns:
            json_output -- directory of trimmed header and end-line in json format
    """

    json_output = []

    for folder in genecluster_json:
        path_output = "/".join(folder.split("/")[:-1]) +"/"
        trim_file = open(path_output+ "genecluster_trim.js", "w")
        json_output += [path_output+"genecluster_trim.js"]
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

    return json_output

def extract_location(json_output):
    """ check the location  of each sequence and keep in dictionary

        Keyword arguments:
            json_output -- directory of trimmed header and end-line in json format (genecluster_trim.js.)

        Returns:
            locus_pHMM -- gene name as key and gene location as value in locus_pHMM dictionary.

            locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary.
     """
    locus_sequence = {}
    locus_location = {}
    locus_pHMM = {}
    for index, dir in enumerate(json_output):
        json_trimfile = json.loads(open(dir).read())
        all_cluster = json_trimfile.keys() # number of cluster in sequences
        # extract all locus tag with location(start,stop)
        for cluster in all_cluster:
            cluster_number = str(cluster)
            for gene in json_trimfile[cluster]['orfs']:
                locus_tag = gene['locus_tag']
                gene = gene['description']
                gene = str(gene)
                location_start = gene.find("Location")
                location_end = gene[location_start:].find("<br>")+location_start
                locus = gene[location_start:location_end].split()
                locus_location[str(locus_tag)] = [locus[1], locus[3]]
        # extract locus tag with signature pHMM
            for gene in json_trimfile[cluster]['orfs']:
                if "Signature pHMM" in gene['description']:
                    locustag = gene['locus_tag']
                    locus_pHMM[str(locustag)] = locus_location[str(locustag)]
                    domain = gene['domains']
                    domain = [dom.replace("plants/", "") for dom in domain]
                    name_domain = ",".join(map(str, domain))
                    details = str(gene['description']).replace("\n", "")
                    total_sequence = re.match("(.*)(copyToClipboard\(')(?P<seq>[A-Z]+)('\))(.*)", details)
                    seq = str(total_sequence.group("seq"))
                    locus_sequence[str(locustag+"|"+name_domain+"|"+cluster_number)] = seq
    # print locus_pHMM
    return locus_pHMM, locus_sequence

def check_overlap(locus_pHMM, locus_sequence):

    """ check the location of each sequence: if the sequence are overlapped, keep the longest sequence.
                                             if the sequence is isoform with the same length, keep the second sequence.

        Keyword arguments:
            locus_pHMM -- gene name  as key and gene location as value in locus_pHMM dictionary.

            locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary

        Returns:
            locus_pHMM  -- gene name as key and gene location as value in locus_pHMM dictionary, with delete overlap location.

            locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary, with delete overlap location.

    """

    for locus_a in locus_pHMM.keys():
        location_a = locus_pHMM[locus_a]
        for locus_b in locus_pHMM.keys():
            location_b = locus_pHMM[locus_b]
            if locus_a == locus_b: # same locus: continue to the next steps
                continue
            else:
                try:
                    if location_a[0] == location_b[0] and location_a[1] == location_b[1]: # comparing start and end location
                        del locus_pHMM[locus_a] # delete first sequence
                    else:
                        if (location_a[0] < location_b[0] < location_a[1]) or (location_a[0] < location_b[1] < location_a[0]):
                            # comparing overlap in case that isoforms do not start at the same position
                            length_1 =  location_a[1] - location_a[0]
                            length_2 =  location_b[1] - location_b[0]

                            if length_1 >= length_2:
                                del locus_pHMM[locus_b]
                            else:
                                del locus_pHMM[locus_a]
                except KeyError:
                        pass

    for locus, sequence in locus_sequence.items():
        if locus.split("|")[0] not in locus_pHMM.keys():
            del locus_sequence[locus]

    return locus_pHMM, locus_sequence



def delseq_isoform(locus_pHMM,locus_sequence, detected_domain, multi_fasta):
    """ Delete the sequence of overlap location, and count the number of genes in cluster_1 and cluster_2

        Keyword arguments:
           locus_pHMM -- gene name as key and gene location as value in locus_pHMM dictionary, with delete overlap location.

           locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary.

           multi_fasta -- directory of multifasta file with delete isoform sequences.

           if -known_domain file.hmm were included in the arguement; file.alined will be written in multi_fasta directory.

        Returns:
            multi_fasta file in multi_fasta directory
            multi_fasta file.aligned.fa in multi_fasta directory

    """

    reversed_domain = ",".join(detected_domain.split(",")[::-1])
    path_to_hmm = hmm_aligned
    for path_domain in path_to_hmm:
        if "\xc2\xa0" in path_domain:
            path_del = path_domain.replace("\xc2\xa0", "") # replace  breaking utf-8 from arguments.
            path_to_hmm +=[path_del]
            path_to_hmm.remove(path_domain)


    output_file = open(multi_fasta, "w")
    for locus_key in locus_sequence.keys():
        # try:
        if locus_key.split("|")[0] not in locus_pHMM.keys():
            del locus_sequence[locus_key]

#write sequence of specific domain without isoform
    locus_keys = locus_sequence.keys()
    locus_keys.sort()
    for locus_key in locus_keys:
        if (str(detected_domain) in locus_key.split("|")[1] or reversed_domain in locus_key.split("|")[1]):
            output_file.write(">{}\n{}\n".format(locus_key, locus_sequence[locus_key]))

    seq_domain = multi_fasta
    output_file.close()

# alignement the extract sequenece with known-domain
    if hmm_aligned == []:
        pass
    else:
        main(seq_domain, path_to_hmm)



if __name__ == "__main__":

    # Defining the fileName variables
    global arguments
    arguments = get_arguments()
    genecluster_json= arguments.input
    multi_fasta = arguments.fasta_output
    detected_domain = arguments.domain
    hmm_aligned = arguments.known_domain

    # Trim and check isoform sequences
    json_output = json_trimheader(genecluster_json)

    #Check overlap location and delete isoform sequences
    locus_pHMM, locus_sequence = extract_location(json_output)
    locus_pHMM, locus_sequence = check_overlap(locus_pHMM, locus_sequence )

    #Delete isoform sequences and write the sequence that contained detected domain
    # if known_domain were included in the argument. The sequence will be aligned with known_domain

    delseq_isoform(locus_pHMM, locus_sequence, detected_domain, multi_fasta)
