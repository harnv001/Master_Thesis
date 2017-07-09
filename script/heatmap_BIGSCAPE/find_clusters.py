#!/usr/bin/env python

"""
Program that match name on phylogenetic tree to the dissimilarity value output from Bigscape program

"""

from __future__ import division

__author__ = "Yosapol Harnvanichvech (950416798110)"
__date__   = "09/06/2016"
__email__  = "yosapol.harnvanichvech@wur.nl"

#Imports
import json
import argparse
import glob


#Functions
def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-name_sequence',
                        dest = 'name_sequence',
                        required = True,
                        type = str,
                        nargs = 1,
                        help = "put name_sequence' in the textfile with one name per line ex. LOCUS|NAME-plant name")
    parser.add_argument('-json_directory',
                        dest ='json_directory',
                        required = True,
                        nargs = '+',
                        help = "put directoy of json file for extract domain information, based on sequence name")
    parser.add_argument('-bigscape_data',
                        dest='bigscape_data',
                        required=True,
                        type=str,
                        nargs=1,
                        help="put the file.network from bigscpae output")
    parser.add_argument('-output_name',
                        dest='output_name',
                        required = True,
                        type = str,
                        nargs = 1,
                        help = "put the file name of output with directory")
    parser.add_argument('-gbk_directory',
                        dest = 'gbk_directory',
                        required = True,
                        type = str,
                        nargs=1,
                        help ="put the heatmap_name")
    parser.add_argument('-seq_selection',
                        dest='seq_selection',
                        required=True,
                        type=str,
                        nargs=1,
                        help="put the name of sequence selection in text file")
    parser.add_argument('-phylo_file',
                        dest='phylo_file',
                        required=True,
                        type=str,
                        nargs=1,
                        help="put the name of phylogenetic tree in newick format")
    #heatmap_name

    #newick file

    return parser.parse_args()

def sequence_name():
    """ parsing full name label from phylogenetic tree and keep in : 1) dictionary of plant name as key and full label as values
                                                                     2) list of only locus name
                                                                     3) list of full label

        Keyword arguments:
            sequence_name -- name of sequence on phylogenetic tree with full label in textfile

        Returns:
            plant_fullname -- dictionary of locu sname as keys and fullname label as values
            locus_name --  list of locus name
            locus_fullname -- list of fullname label
    """
    plant_fullname = {}
    locus_name = []
    locus_fullname = []
    key_name= open(seq_selection[0], "r")
    # print key_name
    for full_name in key_name:
        locus_fullname += [full_name]
        plant_name = full_name.split("|")[0].strip('\n') # tag pattern ex. LOCUS|NAME-plant name
        locus_name += [plant_name]
        plant_fullname[plant_name] = full_name.strip("\n")

    return plant_fullname, locus_name, locus_fullname

def json_trimheader(genecluster_json):
    """Trim header and end line of json file from plantiSMASH output (genecluster.js).

        Keyword arguments:
            genecluster_json -- directory of genecluster.js

        Returns:
            json_output -- directory of trimmed header and end-line in json format
    """

    json_folder = []
    for folder in glob.glob(str(genecluster_json[0])+"/*/genecluster_trim.js"):
        # print dir
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

def locus_cluster(json_folder):
    """ matching cluster order as a key to multiple values of dictionary of locus as a key and domain as a value
        ex. [cluster-009:{AT4G004:TERPENE}

        Keyword arguments:
            json_folder -- directory of trimmed header and end-line in json format (genecluster_trim.js.)

        Returns:
            domain_data -- cluster order as key with multiple values of dictionary of locus as a key and domain as a value
    """
    domain_data = []
    for folder in json_folder:
        json_trimfile = json.loads(open(folder).read())
        all_cluster = json_trimfile.keys() # number of cluster in sequences
        for cluster in all_cluster:
            domains_tmp = {}
            domain_tmp = {}
            for gene in json_trimfile[cluster]['orfs']:
                if "Signature pHMM" in gene['description']:
                    locus_tag = gene['locus_tag']
                    domain = gene['domains']
                    domain = [domain_name.replace("plants/", "") for domain_name in domain] # trim "plant/"
                    name_domain = ",".join(map(str,domain))
                    domain_tmp[str(locus_tag)] = name_domain
            domains_tmp[str(cluster)] = domain_tmp
            domain_data.append(domains_tmp)
    return domain_data

def match_cluster(plant_fullname, locus_name, domain_data):

    """
    Keyword arguments:
        plant_fullname -- sequence_name in the textfile with one name per line"
        locus_name -- list of locus name
        domain_data -- cluster order as key with multiple values of dictionary of locus as a key and domain as a value

    Returns:
        cluster_inf -- list of cluster information which contained : cluster order, plant name and plant fullname(name in phylogeny tree)

    """
    cluster_inf = []
    for all_cluster in domain_data:
        for cluster in all_cluster:
            for locus in all_cluster[cluster]:
                if locus in locus_name:
                    cluster_order = cluster.split("-")
                    cluster_order = "cluster"+str(cluster_order[1].zfill[3]) # make the same format of Genbank file name ex. cluster-009
                    plant_name = plant_fullname[locus].split("-")[-1]
                    cluster_inf += [[cluster_order, plant_name, plant_fullname[locus]]]
    return cluster_inf

def match_locus_gbk(cluster_inf):

    """ matching locus name to Genbank file name (output from plantiSMASH)

        Keyword arguments:
            cluster_inf -- list of cluster information which contained : cluster order, plant name and plant fullname(name in phylogeny tree)

        Returns:
            locus_gbk -- dictionary of Genbank name as a key and locus as values

    """

    locus_gbk = {}
    gbk_folder = glob.glob(str(gbk_directory[0])+"/*/*.gbk") # gbk_directory
    for filename in gbk_folder:
        cluster = filename.split(".")[-2]
        gbk_name = filename.split("/")[-1]
        plant_gbk = filename.split("/")[-2]
        gbk_name = (gbk_name.split(".gbk")[:-1])[0]
        for cluster_detail in cluster_inf:

            if cluster_detail[0] == cluster and cluster_detail[1] == plant_gbk:
                if gbk_name not in locus_gbk.keys():
                    locus_gbk[gbk_name] = [cluster_detail[2]]
                else:
                    locus_gbk[gbk_name].append(cluster_detail[2])
        return locus_gbk

def bigscape_data(bigscape_data):

    """ Read bigscape output.network, delete the header, and keep other information in list

        Keyword arguments:
            bigscape_data -- bigscape outputfile (output.network)

        Returns:
            bigscape_data -- bigscape data in list

    """
    bigscape_data = []
    bigscape_file = open(bigscape_data[0],"r")
    for data in bigscape_file:
        data= data.strip("\n").split("\t")
        bigscape_data += [data]
    del bigscape_data[0] # del header line
    return bigscape_data


def match_data(bigscape_data, locus_gbk):

    """ match cluster distance value from bigscape data to Genbank file name, and change to phylogeny lebel
        Keyword arguments:
            bigscape_data -- bigscape outputfile (output.network)
            locus_gbk -- dictionary of Genbank name as a key and locus as values

        Returns:
            cluster_distance -- dictionary of Genbank file name as keys (tuple of cluster1 and cluster 2) and phylogeny label as values
    """

    cluster_distance = {}
    for index, data in enumerate(bigscape_data):
        cluster_1 = locus_gbk[bigscape_data[index][0]] # cluster 1 = plant full name
        cluster_2 = locus_gbk[bigscape_data][index][1] # cluster 2 = plant full name
        for number_cluster1 in cluster_1:
            for number_cluster2 in cluster_2:
                cluster_distance[number_cluster1, number_cluster2] = locus_gbk[bigscape_data[index][2]] # change Genbank to locus name
                cluster_distance[number_cluster2, number_cluster1] = locus_gbk[bigscape_data[index][2]] # change Genbank to locus name
        return cluster_distance


def seq_selection_distance(cluster_distance,locus_gbk,locus_fullname):

    """ match distance values from Bigscape data to each phylogeny label

            Keyword arguments:
                cluster_distance -- bigscape outputfile (output.network)
                locus_gbk -- dictionary of Genbank name as a key and locus as value
                locus_fullname -- list of fullname label

            Returns:
                distance_matrix -- distance_matrix in the same order as name row and name column

    """
    # seq_select = []
    distance_matrix = []
    # seq_selection = open(seq_selection[0],"r")
    # for seq_name in seq_selection:  # read file and keep name of selected sequence in list
    #     seq_name = seq_name.strip("\n")
    #     seq_select += [seq_name]
    distance_matrix = []
    for name_row in locus_fullname:
        distance_score = []
        for name_col in locus_fullname:
            key_name = name_row, name_col
            if key_name in cluster_distance.keys():
                distance_value = cluster_distance[key_name]
                distance_score +=[str(distance_value)]
            else:
                if [name_col,name_row] or [name_row,name_col] in locus_gbk.values() == True:
                    if name_col == name_row:
                        distance_value = 0 # bigscape did not calculate the values from the same cluster comparison
                        distance_score += [str(distance_value)]
                    else:
                        distance_value = "X"
            distance_matrix += [distance_value]
    return distance_matrix

def write_distance_matrix(distance_matrix, locus_fullname):
    """ write row, column and distance value based on heatmap format from http://itol.embl.de/help.cgi#heatmap

        Keyword arguments:
            distance_matrix -- bigscape outputfile (output.network)
            locus_fullname -- list of fullname label

            Returns:
                  name output.txt for opening in  http://itol.embl.de/help.cgi#heatmap

    """
    matrix_output = open(output_name[0],"w")
    name_col = "\t".join(locus_fullname)
    matrix_output.write("\t" + name_col + "\n")
    for index, name_row in enumerate(locus_fullname):
        distance_score = distance_matrix[index]
        distance_score = "\t".join(distance_score)
        matrix_output.write(name_row + "\t" + distance_score + "\n")





if __name__ == "__main__":
    #Defining the fileName variables
    global arguments
    arguments = get_arguments()
    genecluster_json = arguments.json_directory
    output_name = arguments.output_name #
    gbk_directory = arguments.gbk_directory
    bigscape_data = arguments.bigscape_data
    seq_selection = arguments.name_sequence
    output_name = arguments.output_name
    plant_fullname, locus_name, locus_fullname = sequence_name(seq_selection)
    json_folder = json_trimheader(genecluster_json)
    domain_data = locus_cluster(json_folder)
    cluster_inf = match_cluster(plant_fullname, locus_name, domain_data)
    locus_gbk = match_locus_gbk(cluster_inf)
    bigscape_data = bigscape_data(bigscape_data)
    cluster_distance = match_data(bigscape_data, locus_gbk)
    distance_matrix = seq_selection_distance(cluster_distance, locus_gbk, locus_fullname)
    write_distance_matrix(distance_matrix, locus_fullname)







