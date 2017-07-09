import os
import json
import subprocess
import logging
import sys
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

def main(multi_fasta, path_to_hmm):
    dir_to_fastafile = os.path.join(os.path.dirname(os.path.realpath(__file__)), multi_fasta)
    dir_to_resfile = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                  ("%s%s" % (".".join(multi_fasta.split(".")[:-1]), ".aligned_hmm.fa")))
    dirs_to_hmmfile = [os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                    hmm_file) for hmm_file in path_to_hmm]
    temp_fastas = ""
    results = {}
    with open(dir_to_fastafile, "r") as tpsfa:
        cur_gene = None
        cur_seq = ""
        for curline in tpsfa.readlines():
            curline = curline.rstrip()
            if curline.startswith(">"):
                if cur_gene:
                    temp_fastas += (">%s\n%s\n" % (cur_gene, cur_seq))
                    results[cur_gene] = ""
                cur_gene = "%s" % (curline[1:])
                cur_seq = ""
            else:
                cur_seq += curline
        if cur_gene:
            temp_fastas += (">%s\n%s\n" % (cur_gene, cur_seq))
            results[cur_gene] = ""

    for dir_to_hmmfile in dirs_to_hmmfile:
        hsps = {}
        hit_gids = []
        hmm_len = 0
        with open(dir_to_hmmfile, "r") as hmmfile:
            for line in hmmfile.readlines():
                if line.split(" ")[0] == "LENG":
                    hmm_len = int(line.split(" ")[-1].rstrip())
                    break
        if hmm_len < 1:
            break
        res = run_hmmsearch(dir_to_hmmfile, temp_fastas, 1)
        for runresult in res:
            for hsp in runresult.hsps:
                if hsp.hit.id not in hsps:
                    hsps[hsp.hit.id] = hsp
                if hsps[hsp.hit.id].bitscore < hsp.bitscore:
                    hsps[hsp.hit.id] = hsp
        for gid in hsps:
            hsp = hsps[gid]
            padding_left = ""
            padding_right = ""
            for i in xrange(0, hsp.query_start):
                padding_left += "-"
            for i in xrange(hsp.query_end, hmm_len + 1):
                padding_right += "-"
            if hsp.hit_strand != hsp.query_strand:
                padding_left, padding_right = padding_right, padding_left
            seq = ""
            for c in str(hsp.hit.seq):
                if not c.islower():
                    seq += c
            seq = "".join([padding_left, seq, padding_right])
            results[gid] += seq
            hit_gids.append(gid)
        seq_empty = ""
        for i in xrange(0, hmm_len + 1):
            seq_empty += "-"
        for gid in results:
            if gid not in hit_gids:
                results[gid] += seq_empty

    with open(dir_to_resfile, "w") as resfile:
        for gid in results:
            resfile.write(">%s\n%s\n" % (gid, results[gid]))

def run_hmmsearch(query_hmmfile, target_sequence, cutoff = None):
    "Run hmmsearch"
    command = ["hmmsearch", "--cpu", "2",
               query_hmmfile, '-']
    if cutoff != None:
        if cutoff < 0:
            tc_exist = False
            ga_exist = False
            nc_exist = False
            for line in open(query_hmmfile,"r"):
                firstWord = line.split()[0]
                if firstWord == "TC":
                    tc_exist = True
                if firstWord == "GA":
                    ga_exist = True
                if firstWord == "NC":
                    nc_exist = True
            if tc_exist:
                command.insert(1, "--cut_tc")
            elif ga_exist:
                command.insert(1, "--cut_ga")
            elif nc_exist:
                command.insert(1, "--cut_nc")
            else:
                cutoff = 20
        if cutoff > 0:
            command.insert(1, str(cutoff))
            command.insert(1, "-T")
            command.insert(1, str(cutoff))
            command.insert(1, "--domT")
    out, err, retcode = execute(command, input=target_sequence)
    if retcode != 0:
        print('hmmsearch returned %d: %r while searching %r', retcode,
                        err, query_hmmfile)
        print "error2"
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
        raise

if __name__ == "__main__":
    main()
