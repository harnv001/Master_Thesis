#!/usr/bin/env python

"""
Program that extract the sequence of gene cluster from two input file and
comparing the similarity using heatmap visualization

Requirement: 1. R program
             2. multiple sequence alignemnt MAFFT (MAFFT v7.222 were tested)
"""
from __future__ import division

__author__ = "Yosapol Harnvanichvech (950416798110)"
__date__   = "09/06/2016"
__email__  = "yosapol.harnvanichvech@wur.nl"


#Imports
import json
import re
import os
import subprocess
import argparse
from database import database


#Functions
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-input',
                        dest = 'input',
                        type = str,
                        required = True,
                        nargs = 2,
                        help = "put two genecluster.js file for comparsion")
    parser.add_argument('-cluster',
                        dest='cluster',
                        type=str,
                        nargs= 2,
                        help="cluster number of extraction ex. -cluster cluster-1 cluster-2")
    parser.add_argument('-fasta_output',
                        dest = 'fasta_output',
                        type = str,
                        required = True,
                        help = "put directory of multifasta")
    parser.add_argument('-heatmap_output',
                        dest = 'heatmap_output',
                        type = str,
                        required = True,
                        help = "put directory of heatmap")

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

def extract_location(json_output,cluster_detection):

    """ check the location  of each sequence and keep in dictionary

         Keyword arguments:
             json_output -- directory of genecluster_trim.js.

             cluster_detection -- cluster number for extraction  ex. -cluster cluster-10 cluster-2.

         Returns:
             locus_pHMM -- gene name as key and gene location as value in locus_pHMM dictionary.

             locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary.
      """
    locus_sequence = {}
    locus_location = {}
    locus_pHMM = {}
    # cluster_detection =  arguments.cluster
    for index, dir in enumerate(json_output):
        cluster_number = str(cluster_detection[index])
        json_trimfile = json.loads(open(dir).read())
        all_cluster = json_trimfile.keys() # number of cluster in sequences
        # extract all locus tag with location(start,stop)
        if cluster_number in all_cluster:
            for gene in json_trimfile[cluster_number]['orfs']:
                locus_tag = gene['locus_tag']
                gene = gene['description']
                gene = str(gene)
                location_start = gene.find("Location")
                location_end = gene[location_start:].find("<br>")+location_start
                locus = gene[location_start:location_end].split()
                locus_location[str(locus_tag)] = [locus[1], locus[3]]
        # extract locus tag with signature pHMM
            for gene in json_trimfile[cluster_number]['orfs']:
                if "Signature pHMM" in gene['description']:
                    locustag = gene['locus_tag']
                    locus_pHMM[str(locustag)] = locus_location[str(locustag)]
                    domain = gene['domains']
                    # print domain
                    domain = [dom.replace("plants/", "") for dom in domain]
                    name_domain = ",".join(map(str, domain))
                    details = str(gene['description']).replace("\n", "")
                    total_sequence = re.match("(.*)(copyToClipboard\(')(?P<seq>[A-Z]+)('\))(.*)", details)
                    seq = str(total_sequence.group("seq"))
                    locus_sequence[str(locustag+"|"+name_domain+"|"+cluster_number)] = seq
    # delete genecluster_trim.js
    for folder in json_output:
        cmd = "rm {}".format(folder)
        if subprocess.check_call(cmd, shell=True) == 0:
            print ("delete"+ str(folder))
    return locus_pHMM, locus_sequence

# delete the isoforms by checking overlapped location
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

    for locus, sequence in locus_sequence.items():
        if locus.split("|")[0] not in locus_pHMM.keys():
            del locus_sequence[locus]

    return locus_pHMM, locus_sequence

#delete isoform key
def delseq_isoform(locus_sequence, cluster_detection, multi_fasta):

    """ Delete the sequence of overlap location, and count the number of genes in cluster_1 and cluster_2

       Keyword arguments:
          locus_sequence -- gene name as key and amino acid sequence as value in locus_pHMM dictionary.

          cluster_detection -- the order of cluster that select for detection

          multi_fasta -- directory of multifasta file with delete isoform sequences.

       Returns:
           count_gene_cluster1 -- number of gene in cluster 1.

           count_gene_cluster2  --  number of gene in cluster 2.

    """
    number_cluster1 = cluster_detection[0]
    number_cluster2 = cluster_detection[1]
    count_gene_cluster1 = 0
    count_gene_cluster2 = 0
    output_file = open(multi_fasta, "w")
    locus_key = locus_sequence.keys()
    locus_key.sort()
    for locus in locus_key:
        if locus.split ("|")[-1] == number_cluster1:
            count_gene_cluster1 +=1
        if locus.split("|")[-1] == number_cluster2:
            count_gene_cluster2 +=1

        if locus.split("|")[1] in database.keys():
            domain = database[locus.split("|")[1]][0]
            label_name = str(locus.split("|")[0]+"|"+ domain)
            print (label_name)
        else:
            label_name = locus.split("|cluster")[0]
            print (label_name)

        output_file.write(">{}\n{}\n".format(label_name, locus_sequence[locus]))
    output_file.close()
    return count_gene_cluster1, count_gene_cluster2

# run alignment using MAFFT v7.222.
def run_alignment(multi_fasta):

    """ Runs multiple sequence alignment program (MAFFT v7.222) on the command line using subprocess

       Keyword arguments:
           multi_fasta -- directory of multifasta file with delete isoform sequences.

       Returns:
           fasta_aligned -- output directory of sequence alignment.

   """
    fasta_aligned = multi_fasta.split(".")[:-1][0] + "_aligned.fa"
    cmd = ["mafft.bat","--auto",multi_fasta]
    if not os.path.isfile(fasta_aligned):
        with open(fasta_aligned,'w') as out_file:
            process = subprocess.Popen(cmd, stdout=out_file)
            process.communicate()
            return_code = process.returncode
            if return_code == 0:
                print ("run alignement from MAFFT v7.222")
            else:
                print ("run alignement error")
    else:
        print (fasta_aligned.split("/")[-1] + " already exists!")
    return fasta_aligned

def heatmap_similarity(fasta_aligned , count_gene_cluster1, count_gene_cluster2, heatmap_output):

    """ Runs heatmap similarity from seqINR package (R program) on the command line using subprocess

        Keyword arguments:
            fasta_aligned -- output directory of sequence alignment

            count_gene_cluster1 -- number of gene in cluster 1

            count_gene_cluster2  --  number of gene in cluster 2

            heatmap_output -- directory of heatmap figure in .png file

        Returns:
            fasta_aligned -- output directory of sequence alignedment

    """

    total_cluster = count_gene_cluster1 + count_gene_cluster2
    cmd = ['Rscript', 'cal_similarity.R', fasta_aligned, str(count_gene_cluster2), str(total_cluster), heatmap_output]
    if not os.path.isfile(heatmap_output):
        process = subprocess.Popen(cmd)
        process.communicate()
        return_code = process.returncode
        if  return_code == 0:
            print ("finish cluster heatmap comparison")
        else:
            print ("run heatmap error")
    else:
        print (heatmap_output.split("/")[-1] + " already exists!")

if __name__ == "__main__":

    #Defining the fileName variables
    global arguments
    arguments = get_arguments()
    genecluster_json= arguments.input
    multi_fasta = arguments.fasta_output
    cluster_detection = arguments.cluster
    heatmap_output = arguments.heatmap_output

    #Trim and check isoform sequences
    json_output = json_trimheader(genecluster_json)
    locus_pHMM, locus_sequence = extract_location(json_output,cluster_detection)

    #Check overlap location and delete isoform sequences
    locus_pHMM,locus_sequence = check_overlap(locus_pHMM,locus_sequence)
    count_gene_cluster1, count_gene_cluster2 = delseq_isoform(locus_sequence, cluster_detection,multi_fasta)

    #Run multiple sequence alignemnt using MAFFT v7.222
    fasta_aligned = run_alignment(multi_fasta)

    #Comparing heatmap similarity between genecluster using seqINR package in R program
    heatmap_similarity(fasta_aligned, count_gene_cluster1, count_gene_cluster2,heatmap_output)