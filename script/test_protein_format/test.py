import json
import re
import os
import glob
import sys
from Bio import SeqIO


# match protein tag from Bigscape to phylogenetic name

def jstrim_directory():
    jstrim_dir = []
    for genecluster_json in glob.glob("/home/harnv001/public_html/TPS_delisoform/*/geneclusters_trim.js"):
        jstrim_dir += [genecluster_json]
    return jstrim_dir


#
# clade name from phylogenetic: matching shortname with longname
def sequence_selection(file):
    dict_fullname = {}
    list_key = []
    key_name = open(file, "r")
    # print key_name
    for full_name in key_name:
        # print full_name
        name = full_name.split("|")[0].strip('\n')
        list_key += [name]
        dict_fullname[name] = full_name.strip("\n")

    return dict_fullname, list_key


# get domain in each cluster
def check_cluster(jstrim_dir, list_key):
    plants_domains = []
    for dir in jstrim_dir:
        # print dir
        json_trimfile = json.loads(open(str(dir)).read())
        # print json_trimfile
        total_cluster = json_trimfile.keys()
        # print total_cluster
        database = {"p450": ["Cytochrome_450", "#ff0000", "RE"],
                    "Terpene_synth_C": ["Terpene synthase", "#35F800", "RE"],
                    "SQHop_cyclase_C": ["Terpene synthase", "#35F800", "RE"],
                    "SQHop_cyclase_N": ["Terpene synthase", "#35F800", "RE"],
                    "SQHop_cyclase_C,SQHop_cyclase_N": ["Terpene synthase", "#35F800", "RE"],
                    "Terpene_synth": ["Terpene synthase", "#35F800", "RE"],
                    "Lycopene_cycl": ["Terpene synthase", "#35F800", "RE"],
                    "Cu_amine_oxid": ["Copper amine oxidase", "#E9C63E", "RE"],
                    "Bet_v_1": ["Pictet-Spengler enzyme (Bet v1)", "#735100", "RE"],
                    "Glycos_transf_1": ["Glycosyltransferase", "#F477A6", "RE"],
                    "Glycos_transf_2": ["Glycosyltransferase", "#F477A6", "RE"],
                    "Glyco_transf_28": ["Glycosyltransferase", "#F477A6", "RE"],
                    "UDPGT": ["Glycosyltransferase", "#F477A6", "RE"],
                    "UDPGT_2": ["Glycosyltransferase", "#F477A6", "RE"],
                    "Chal_sti_synt_C": ["Ketosynthase", "#75FFD8", "RE"],
                    "Chal_sti_synt_N": ["Ketosynthase", "#75FFD8", "RE"],
                    "SE": ["Squalene epoxidase", "#206B14", "RE"], "COesterase": ["COesterase", "#A221FF", "RE"],
                    "Methyltransf_2": ["Methyltransferase", "#CF0FFF", "RE"],
                    "Methyltransf_3": ["Methyltransferase", "#CF0FFF", "RE"],
                    "Methyltransf_7": ["Methyltransferase", "#CF0FFF", "RE"],
                    "Methyltransf_11": ["Methyltransferase", "#CF0FFF", "RE"],
                    "cMT": ["Methyltransferase", "#CF0FFF", "RE"], "nMT": ["Methyltransferase", "#CF0FFF", "RE"],
                    "oMT": ["Methyltransferase", "#CF0FFF", "RE"],
                    "Transferase": ["BAHD acyltransferase", "#0003EC", "RE"],
                    "Peptidase_S10": ["Scl acyltransferase", "#0c126d", "RE"],
                    "Epimerase": ["Epimerase", "#F245BD", "RE"],
                    "adh_short": ["Oxidoreductase", "#4DD0E1", "RE"],
                    "adh_short_C2": ["Oxidoreductase", "#4DD0E1", "RE"],
                    "NAD_binding_1": ["Oxidoreductase", "#4DD0E1", "RE"],
                    "GMC_oxred_N": ["Oxidoreductase", "#4DD0E1", "RE"],
                    "GMC_oxred_C": ["Oxidoreductase", "#4DD0E1 ", "RE"], "DIOX_N": ["Dioxygenase", "#F5B2FF", "RE"],
                    "2OG-FeII_Oxy": ["Dioxygenase", "#F5B2FF", "RE"],
                    "AMP-binding": ["CoA-ligase", "#8F400B", "RE"], "Amino_oxidase": ["Amino oxidase", "#2C113A", "RE"],
                    "Aminotran_1_2": ["Aminotransferase", "#E42E00", "RE"],
                    "Aminotran_3": ["Aminotransferase", "#E42E00", "RE"],
                    "Prenyltrans": ["Prenyltransferase", "#552288", "RE"],
                    "Prenyltransf": ["Prenyltransferase", "#552288", "RE"],
                    "UbiA": ["Prenyltransferase", "#552288", "RE"],
                    "Str_synth": ["Strictosidine synthase-like", "#dd71dd", "RE"],
                    "PRISE": ["PRISE enzymes", "#13C600", "RE"], "Dirigent": ["Dirigent enzymes", "#606000", "RE"]
                    }

        for cluster in total_cluster:
            locus_domains = {}  # each plants
            locus_domain = {}
            name_cluster = str(cluster)
            # print name_cluster
            for gene in json_trimfile[cluster]['orfs']:
                if "Signature pHMM" in gene['description']:
                    locustag = gene['locus_tag']

                    domain = gene['domains']
                    domain = [dom.replace("plants/", "") for dom in domain]
                    name_domain = ",".join(map(str, domain))

                    locus_domain[str(locustag)] = name_domain

            locus_domains[str(cluster)] = locus_domain
            plants_domains.append(locus_domains)
    return plants_domains
    # print plants_domains


# extract cluster, plant name, full name
def match_cluster(list_key, dict_fullname, plants_domains):
    # file = open("/home/harnv001/public_html/find_cluster.txt","w")
    cluster_details = []
    count = 0
    #     #print locus_domains
    for all_cluster in plants_domains:
        for cluster in all_cluster:
            check = False
            for locus_key in all_cluster[cluster]:
                # print locus_key
                if locus_key in list_key:
                    count += 1
                    number_cluster = cluster.split("-")
                    number_cluster = number_cluster[1].zfill(3)
                    plant_name = dict_fullname[locus_key].split("-")[-1]
                    # print plant_name
                    cluster_details += [['cluster' + str(number_cluster), plant_name, dict_fullname[locus_key]]]
                    # print number_cluster
                    # file.write(str(count)+ '\t'+ cluster+'\t'+ dict_fullname[locus_key]+'\t' + locus_key + '\t' + all_cluster[cluster][locus_key]+'\n')
    # print cluster_details
    return cluster_details

    # match TPS name with gbk output file name from plantiSMASH


def match_TPS_gbk(cluster_detials):
    gbk_TPS = {}  # match gbk and TPS name
    gbk_dir = glob.glob("/home/harnv001/public_html/gbk/*/*.gbk")
    for filename in gbk_dir:
        cluster = filename.split(".")[-2]
        gbk_name = filename.split("/")[-1]
        gbk_name = (gbk_name.split(".gbk")[:-1])[0]
        plant_gbk = filename.split("/")[-2]
        for detail in cluster_details:
            if detail[0] == cluster and detail[1] == plant_gbk:
                if gbk_name not in gbk_TPS.keys():
                    gbk_TPS[gbk_name] = [detail[2]]
                else:
                    gbk_TPS[gbk_name].append(detail[2])
    # print gbk_TPS
    # print gbk_TPS
    return gbk_TPS


# extract data from bigscape file in each linek
def extract_bigscape():
    bigscape_data = []
    bigscape_file = open("/home/harnv001/public_html/final_result/non_core_new/non_core_new_all_mix_c1.00.txt", "rU")
    for line in bigscape_file:
        # print line
        line = line.strip('\n').split("\t")
        # print line
        bigscape_data += [line]
    del bigscape_data[0]  # del header line
    return bigscape_data
    # print bigscape_data


def match_data(bigscape_data, gbk_TPS):
    cluster_distance = {}
    # change gbk name to TPS name
    for index, data in enumerate(bigscape_data):
        # print index,data
        cluster_1 = gbk_TPS[bigscape_data[index][0]]  # cluster 1 = plant full name
        cluster_2 = gbk_TPS[bigscape_data[index][1]]  # cluster 2 = plant full name
        for number_cluster1 in cluster_1:
            for number_cluster2 in cluster_2:
                cluster_distance[number_cluster1, number_cluster2] = bigscape_data[index][2]
                cluster_distance[number_cluster2, number_cluster1] = bigscape_data[index][2]
                # #         # break
    return cluster_distance


#     print cluster_distance
def mapdata_matrix(cluster_distance, gbk_TPS):
    # print gbk_TPS
    TPS_name = []
    score_matrix = []
    TPS_file = open("/home/harnv001/public_html/final_result/complate_plantname.txt", "r")
    # TPS_file = open("/home/harnv001/public_html/")
    for TPS_line in TPS_file:
        TPS_line = TPS_line.strip('\n')
        TPS_name += [TPS_line]
    # loop TPS with TPS
    for TPS_row in TPS_name:
        score = []
        for TPS_col in TPS_name:
            key = TPS_row, TPS_col
            # print key
            # print cluster_distance
            if key in cluster_distance.keys():
                distance_value = cluster_distance[key]
                score += [str(distance_value)]
            else:
                # if [TPS_row,TPS_col] or [TPS_col,TPS_row] in gbk_TPS.values() == True: #same cluster were not calculate from Bigscape; check if the plant are the same cluster, put 0:
                if TPS_row == TPS_col or TPS_col == TPS_row:
                    distance_value = 0
                    score += [str(distance_value)]
                else:
                    score += ["1"]
            score_matrix += [score]

        # print score_matrix
        # break
    return score_matrix, TPS_name
    # return score_matrix, TPS_name


#
def write_matrix(score_matrix, TPS_name):
    #     # print TPS_name
    #     matrix = open("/home/harnv001/public_html/gbk/core_biosynthetic/gbkmatrix_TPS_TPS-no-anchore_domain_test.txt", "w") # edit name output
    #    matrix = open("/home/harnv001/public_html/final_result/non_core_merge_all_mix_c1.00/result_non_core_merge_all_mix_c1.00.txt", "w")
    matrix = open("/home/harnv001/public_html/final_result/non_core_new/result_non_core_new2.txt", "w")

    TPS_col = '\t'.join(TPS_name)
    #
  x  #     # print ("\t" + TPS_col +"\n")
    matrix.write("\t" + TPS_col + "\n")
    # c
    for index, TPS_row in enumerate(TPS_name):
        score_row = score_matrix[index]
        score = '\t'.join(score_row)
        matrix.write(TPS_row + "\t" + score + "\n")


if __name__ == "__main__":
    jstrim_dir = jstrim_directory()
    dict_fullname, list_key = sequence_selection(file="/home/harnv001/public_html/TPS_delisoform/clade_selection.txt")
    plants_domains = check_cluster(jstrim_dir, list_key)
    cluster_details = match_cluster(list_key, dict_fullname, plants_domains)
    gbk_TPS = match_TPS_gbk(cluster_details)
    bigscape_data = extract_bigscape()
    cluster_distance = match_data(bigscape_data, gbk_TPS)
    score_matrix, TPS_name = mapdata_matrix(cluster_distance, gbk_TPS)
    write_matrix(score_matrix, TPS_name)
