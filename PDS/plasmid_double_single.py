#Counter-marker-assisted-HR-based_dsDNA HS_point mutation_Gibson assembly

import os
import primer3
import json
import pandas as pd
import configparser
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Obtain the reverse complementary sequence 
def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

def conf_read(filename): 
    config = configparser.ConfigParser()
    config.read(filename)
    res = dict(config._sections["point"])
    return res

def blast_search(input_file_path,genome,workdir):
    """
    Arguments:
        input_file_path[str]:  input information
        genome[str]:  edit reference
        workdir[str]: dir path
    Return [dict]
        {
            "key1":{
                "chrom":"NC_001133.9",
                "start":"1000",
                "unique_mapped":1,  次数 100%比对
                "mapped":1,  次数
                "reverse_mapped":1, 次数
                "description":"NC_001133.9:26529-26708;NC_001133.9:26664-26843;",  
            },

        }

    """
    blast_output_file_path=os.path.join(workdir,'blast_search.txt')

    ref_lib=genome.split('/')[-1].split('.')[0]
    input_fasta = os.path.join(workdir,'blast.fasta')
    fasta_length_dict = {}
    my_records = []
    with open(input_file_path,'r') as ifile:
        index = 0 
        for line in ifile:
            if index != 0 :
                linelist = line.split(',')
                seqid = linelist[0]
                seq = linelist[1].strip()
                rec = SeqRecord(Seq(seq),id=seqid)

                fasta_length_dict[seqid] = len(seq)
                my_records.append(rec)
            index += 1
    # input fasta
    SeqIO.write(my_records, input_fasta, "fasta")
    # run blast
    os.system("makeblastdb -in "+genome+" -dbtype nucl -parse_seqids -out "+ref_lib)
    os.system("blastn -query "+input_fasta+" -db "+ref_lib+" -outfmt 6 -out "+blast_output_file_path+" -evalue 1e-30 ")

    # return
    dictall = {}
    with open(blast_output_file_path,"r") as f:
        for i in f:
            linelist = i.split('\t')
            key = linelist[0]
            chrom = linelist[1]
            identity = linelist[2]
            allength = linelist[3]
            start = linelist[8]
            end = linelist[9]

            if key not in dictall:
                dictall[key] = {
                    "chrom":chrom,
                    "start":start,
                    "unique_mapped":1 if (int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key]) else 0,
                    "mapped":1,
                    "reverse_mapped":1 if (int(start) > int(end) and int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key]) else 0,
                    "description":'%s:%s-%s;' %(chrom,start,end),
                }
            else:
                dictall[key]["mapped"] += 1
                if int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key] :
                    dictall[key]["unique_mapped"] += 1
                    dictall[key]["chrom"] = chrom
                    dictall[key]["start"] = start
                    if int(start) > int(end):
                        dictall[key]["reverse_mapped"] += 1
                    dictall[key]["description"] += '%s:%s-%s;' %(chrom,start,end)
    return dictall

def input_to_primer_template(input_file_path,genome,config,workdir):
    """
    Arguments:
        input_file_path[str]:  input information
        genome[str]:  edit reference
        config[dic]: points information
    Return [dict]
        {
            "key1":{
                "seq_uha_max_whole":"",
                "seq_dha_max_whole":"",
                "seq_altered":"",
                "seq_dr":"",
                "type":"",   # [substitution,deletion,insertion]
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_max_whole  down 100bp sequence,
            },
            "key2":{
                "seq_uha_max_whole":"",
                "seq_dha_max_whole":"",
                "seq_altered":"",
                "seq_dr_uha":"",
                "seq_dr_dha":"",
                "type":"",
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_max_whole  down 100bp sequence,
            }
        }
    """
    max_left_arm_seq_length=int(config['max_left_arm_seq_length'])
    max_right_arm_seq_length=int(config['max_right_arm_seq_length'])
    # 如果未上传质粒，max_verify_1_up_ponit max_verify_2_down_ponit为0
    max_verify_2_up_ponit=int(config['max_verify_2_up_ponit'])
    max_verify_2_down_ponit=int(config['max_verify_2_down_ponit'])
    dr_seq_length=int(config['dr_seq_length'])
    length_threshold_msddsc=int(config['length_threshold_msddsc'])
    input_format=input_file_path.split('.')[-1]
    primer_template = {}
    if input_format != 'csv':
        error_message = "The input file format needs to be 'csv'"
        print(error_message)
        return error_message
    else:
        # read genome
        # record = SeqIO.read(genome, "fasta").seq
        # record_list = list(SeqIO.parse(genome, "fasta"))
        record_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
        blast_search_dict = blast_search(input_file_path,genome,workdir)
        df=pd.read_csv(input_file_path)
        num_lines_df=df.shape[0]
        for i in range(num_lines_df):
            data=df.loc[i].values
            #print(data)
            mun_id=str(data[0])
            if len(data) < 5:
                error_message = "Some necessary input information is missing in the line "+mun_id+" of the input file"
                print(error_message)
                return  error_message
            else:
                upstream = data[1].strip().upper()
                ref = data[2]
                alt = data[3]
                mutation_type = data[4].strip().lower()
                
                if mun_id not in blast_search_dict:
                    # 没有blast上
                    error_message = "The upstream sequence of " + mun_id + " can not be mapped to the target genome. please check whether the target sequence is located on the target genome. Please rightly prepare input file for target manipulation as the example of 2,3-BD."
                    print(error_message)
                    return error_message
                elif blast_search_dict[mun_id]["unique_mapped"] > 1:
                    # 多次比对上 100
                    error_message = "The upstream sequence of " + mun_id + "  can be mapped to multiple loci in the target genome, %s, Please provide a longer upstream seqeunce." % blast_search_dict[mun_id]["description"]
                    print(error_message)
                    return error_message
                elif blast_search_dict[mun_id]["unique_mapped"] == 0:
                    # 无 100 比对
                    error_message = "The upstream sequence of " + mun_id + " can not be uniquely mapped to the target genome. Please check whether the target sequence is located on the target genome."
                    print(error_message)
                    return error_message
                elif blast_search_dict[mun_id]["unique_mapped"] == 1:
                    # 开始突变的碱基在genome上的索引
                    if blast_search_dict[mun_id]["reverse_mapped"]:
                        record = revComp(str(record_dict[blast_search_dict[mun_id]["chrom"]].seq))
                        upstream_start_index = len(record) - int(blast_search_dict[mun_id]["start"])
                        # record=str(record_dict[blast_search_dict[mun_id]["chrom"]].seq).translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
                        # mutation_pos_index = int(blast_search_dict[mun_id]["start"])
                        strand = "minus"
                    else:
                        record = str(record_dict[blast_search_dict[mun_id]["chrom"]].seq)
                        upstream_start_index = int(blast_search_dict[mun_id]["start"])-1
                        strand = "plus"
                    mutation_pos_index = upstream_start_index + len(upstream)
                    mutation_pos_ref = record[mutation_pos_index].upper()
                    # length 
                    if mutation_pos_index - max_left_arm_seq_length - max_verify_2_up_ponit < 0:
                        error_message = "The length of upstream sequence of manipulation site of " + mun_id + " must be larger than sum of 'Max Length of UHA' and 'Max Length of UIS'."
                        print(error_message)
                        return error_message
                    if mutation_type == "deletion":
                        # 先判断用户提供的ref是否与参考基因组ref一致
                        genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]
                        if genome_ref.upper() == ref.upper():
                            primer_template[mun_id] = {
                                    "seq_uha_max_whole":str(record[
                                        mutation_pos_index - max_left_arm_seq_length - max_verify_2_up_ponit          :mutation_pos_index
                                        ]),
                                    "seq_dha_max_whole":str(record[
                                        mutation_pos_index + len(ref)
                                        : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                        ]),
                                    "seq_dr":str(record[
                                        mutation_pos_index - max_left_arm_seq_length:mutation_pos_index
                                        ])[-dr_seq_length:],
                                    "seq_altered": "",
                                    "ref": ref,
                                    "strand":strand,
                                    "mutation_pos_index":mutation_pos_index,
                                    "geneid":blast_search_dict[mun_id]["chrom"],
                                    "uha_upstream":str(record[
                                        mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit -100
                                        : mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit
                                        ]),
                                    "dha_downstream":str(
                                        record[
                                            mutation_pos_index + max_right_arm_seq_length + max_verify_2_down_ponit
                                            : mutation_pos_index + max_right_arm_seq_length + max_verify_2_down_ponit + 100
                                        ]),
                                    "level":"del",
                                    "type":mutation_type,
                                }
                    elif mutation_type in ["insertion","substitution"]:
                        if mutation_type == "insertion":
                            ref = ""
                        genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]
                        if genome_ref.upper() == ref.upper():
                            if len(alt)<=dr_seq_length:
                                primer_template[mun_id] = {
                                    "seq_uha_max_whole":str(record[
                                        mutation_pos_index - max_left_arm_seq_length - max_verify_2_up_ponit:mutation_pos_index
                                        ]),
                                    "seq_dha_max_whole":str(record[
                                        mutation_pos_index + len(ref)
                                        : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                        ]),
                                    "seq_altered":alt,
                                    "ref": ref,
                                    "strand":strand,
                                    "mutation_pos_index":mutation_pos_index,
                                    "geneid":blast_search_dict[mun_id]["chrom"],
                                    "seq_dr":str(record[
                                        mutation_pos_index - max_left_arm_seq_length:mutation_pos_index
                                        ]+alt)[-dr_seq_length:],
                                    "uha_upstream":str(record[
                                        mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit -100
                                        : mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit

                                        ]),
                                    "dha_downstream":str(
                                        record[
                                            mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                            : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit + 100
                                        ]),
                                    "level":"<30",
                                    "type":mutation_type,
                                }
                            elif dr_seq_length<len(alt)<=length_threshold_msddsc:
                                primer_template[mun_id] = {
                                    "seq_uha_max_whole":str(record[
                                        mutation_pos_index - max_left_arm_seq_length - max_verify_2_up_ponit :mutation_pos_index
                                        ]),
                                    "seq_dha_max_whole":str(record[
                                        mutation_pos_index + len(ref)
                                        : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                        ]),
                                    "seq_altered" : alt,
                                    "ref": ref,
                                    "strand":strand,
                                    "mutation_pos_index":mutation_pos_index,
                                    "geneid":blast_search_dict[mun_id]["chrom"],
                                    "seq_dr_uha":alt[:math.ceil((len(alt)+dr_seq_length)/2)],
                                    "seq_dr_dha":alt[math.ceil((len(alt)-dr_seq_length)/2):],
                                    "uha_upstream":str(record[
                                        mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit -100
                                        : mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit
                                        ]),
                                    "dha_downstream":str(
                                        record[
                                            mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                            : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit + 100
                                        ]),
                                    "level":">30&<90",
                                    "type":mutation_type,
                                }
                            else:
                                primer_template[mun_id] = {
                                    "seq_uha_max_whole":str(record[
                                        mutation_pos_index - max_left_arm_seq_length - max_verify_2_up_ponit :mutation_pos_index
                                        ]),
                                    "seq_dha_max_whole":str(record[
                                        mutation_pos_index + len(ref)
                                        : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                        ]),
                                    "seq_altered" : alt,
                                    "ref": ref,
                                    "strand":strand,
                                    "mutation_pos_index":mutation_pos_index,
                                    "geneid":blast_search_dict[mun_id]["chrom"],
                                    "seq_dr":alt[-dr_seq_length:],
                                    "uha_upstream":str(record[
                                        mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit -100
                                        : mutation_pos_index 
                                        - max_left_arm_seq_length - max_verify_2_up_ponit
                                        ]),
                                    "dha_downstream":str(
                                        record[
                                            mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit
                                            : mutation_pos_index + len(ref) + max_right_arm_seq_length + max_verify_2_down_ponit + 100
                                        ]),
                                    "level":">90", # 增加一对引物
                                    "type":mutation_type,
                                }
                    else:
                        error_message = "The target manipulation type of " + mun_id + " must be equal to 'insertion,substitution or deletion', Please rightly prepare input file for target manipulation as the example of 2,3-BD."
                        print(error_message)
                        return  error_message
    return primer_template

#Obtain linear plasmid sequence
#linear plasmid and screening marker seq both use this
def linearSeq(input_file_path):
    with open(input_file_path,'r') as input_file:
        return input_file.read().upper()

#Primer design for upstream Homologous Arm
def primerDesign(seqId,seqTemplate,config,stype):
    """
    Description: primer design
    Args:
        seqId ([str]): [sequence id from input file]
        seqTemplate ([str]): [sequence]
        config ([dict]): [points dict]
        stype ([str]): ["verify_1","verify_2","left_arm","right_arm"]

    Returns:
        [dict]: [output of primer3]
    """
    seqlength=len(seqTemplate)
    seq_args = {
            'SEQUENCE_ID': seqId,
            'SEQUENCE_TEMPLATE': seqTemplate,
        }
    if stype =="verify_2":
        print("************* %s process ************" %stype)
        left_length=int(config["max_verify_2_up_ponit"])-int(config["min_verify_2_up_ponit"])
        right_length=int(config["max_verify_2_down_ponit"])-int(config["min_verify_2_down_ponit"])
        region_list = [[0,left_length,seqlength-right_length-1,right_length]]
        size_range = [seqlength-left_length-right_length+36,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    elif stype == "left_arm":
        print("************* %s process ************" %stype)
        left_length=int(config["max_left_arm_seq_length"])-int(config["min_left_arm_seq_length"])
        size_range = [seqlength-left_length+18,seqlength]
        region_list = [[0,left_length+18,seqlength-25,25]]
        seq_args["SEQUENCE_FORCE_RIGHT_START"] = seqlength-1
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    elif stype == "right_arm":
        print("************* %s process ************" %stype)
        right_length=int(config["max_right_arm_seq_length"])-int(config["min_right_arm_seq_length"])
        size_range = [seqlength-right_length+18,seqlength]
        region_list = [[0,25,seqlength-right_length-1,right_length]]
        seq_args["SEQUENCE_FORCE_LEFT_START"] = 0
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    elif stype == "insert" or "screening_marker":
        print("************* %s process ************" %stype)
        size_range = [seqlength,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = [[0,25,seqlength-26,25]]
        seq_args["SEQUENCE_FORCE_LEFT_START"] = 0
        seq_args["SEQUENCE_FORCE_RIGHT_START"] = seqlength-1
    elif stype == "verify_1":
        print("************* %s process ************" %stype)
        left_length=int(config["max_verify_1_up_ponit"])-int(config["min_verify_1_up_ponit"])
        right_length=int(config["max_verify_1_down_ponit"])-int(config["min_verify_1_down_ponit"])
        region_list = [[0,left_length,seqlength-right_length-1,right_length]]
        size_range = [seqlength-left_length-right_length+36,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    global_args = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': float(config[stype+"_primer_opt_tm"]),
            'PRIMER_MIN_TM': float(config[stype+"_primer_min_tm"]),
            'PRIMER_MAX_TM': float(config[stype+"_primer_max_tm"]),
            'PRIMER_MIN_GC': float(config[stype+"_primer_min_gc"]),
            'PRIMER_MAX_GC': float(config[stype+"_primer_max_gc"]),
            'PRIMER_PICK_ANYWAY':1,
            'PRIMER_PRODUCT_SIZE_RANGE': size_range,
            'PRIMER_NUM_RETURN':1
    }
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
    return primer3_result

#Fill in the "successful" and "uha" elements of dict_primers_whole based on the results of the upstream homologous sequence primer design to obtain the downstream homologous sequence primer design template and seq_uha_left_point
def leftPrimerAttribute(ID,dict_left_primers,target_seq,dict_input_seq,screeningmarker_seq,config,plasmidseq):
    dict_left_primers_attribute={}
    length_homo_seq=int(config['length_homologous_sequence'])
    max_verify_2_up_ponit=int(config['max_verify_2_up_ponit'])
    dict_left_primers_attribute['ID']=ID
    uha_left_primer=dict_left_primers['PRIMER_LEFT_0_SEQUENCE']
    uha_right_primer=dict_left_primers['PRIMER_RIGHT_0_SEQUENCE']
    uha_right_primer_rev=revComp(uha_right_primer)
    uha_temp=dict_input_seq['seq_uha_max_whole'][len(dict_input_seq['seq_uha_max_whole'])-int(config['max_left_arm_seq_length']):]
    seq_uha_left_point_on_temp=uha_temp.find(uha_left_primer)
    # 130 增加序列长度，提高stringfind的精度
    seq_uha_left_primer_add_130=uha_temp[seq_uha_left_point_on_temp:seq_uha_left_point_on_temp+130]
    seq_uha_left_point=target_seq.find(seq_uha_left_primer_add_130)
    seq_uha_right_point=uha_temp.find(uha_right_primer_rev)
    #uha_left_primer_add_seq=plasmidseq[0-length_homo_seq:]
    # P5前20bp
    uha_right_primer_add_seq=screeningmarker_seq[:length_homo_seq]
    uha_left_primer_whole = plasmidseq[-length_homo_seq:] + uha_left_primer
    dict_left_primers_attribute['PRIMER_LEFT_WHOLE_SEQUENCE']=uha_left_primer_whole
    if  dict_input_seq["level"] == ">30&<90":
        # 处理insertion substitution的 30< alt <90
        rev_uha_right_primer=uha_right_primer_rev+dict_input_seq["seq_dr_uha"]+uha_right_primer_add_seq
    elif dict_input_seq["level"] == ">90":
        # 处理insertion substitution的 alt >90
        rev_uha_right_primer=uha_right_primer_rev+ dict_input_seq["seq_altered"][:length_homo_seq]
    else:
        # 处理deletion  和 insertion substitution的 alt <= 30
        rev_uha_right_primer=uha_right_primer_rev+dict_input_seq["seq_altered"]+uha_right_primer_add_seq
    seq_uha=plasmidseq[-length_homo_seq:] + uha_temp[seq_uha_left_point_on_temp:seq_uha_right_point]+rev_uha_right_primer
    dict_left_primers_attribute['PRIMER_RIGHT_WHOLE_SEQUENCE']=revComp(rev_uha_right_primer)
    seq_uha_without_add=uha_temp[seq_uha_left_point_on_temp:seq_uha_right_point]+uha_right_primer_rev
    dict_left_primers_attribute['SEQ_UHA_WITHOIUTE_ADD']=seq_uha_without_add
    seq_uha_left_point_on_seq_uha_max_whole=dict_input_seq['seq_uha_max_whole'].find(uha_left_primer)
    seq_verify2_temp_uha_half= dict_input_seq['seq_uha_max_whole'][seq_uha_left_point_on_seq_uha_max_whole-max_verify_2_up_ponit:seq_uha_left_point_on_seq_uha_max_whole]+seq_uha_without_add
    dict_left_primers_attribute['SEQ_VERIFY2_TEMP_UHA_HALF']=seq_verify2_temp_uha_half
    #dha_temp_add=dict_input_seq['seq_uha_max_whole'][-50:]+dict_input_seq['seq_altered']+dict_input_seq['seq_dha_max_whole'][:int(config['max_right_arm_seq_length'])]
    #dha_5end=dha_temp_add.find(rev_uha_right_primer)+len_rev_uha_right_primer-1
    #dha_temp=dha_temp_add[dha_5end+1:]
    dict_left_primers_attribute['PRODUCT_SEQUENCE']=seq_uha
    #type3_left_primer_add_seq=seq_uha[0-int(config['length_homologous_sequence']):]
    dict_left_primers_attribute['PRODUCT_WHOLE_LENGTH']=len(seq_uha)
    return dict_left_primers_attribute
    #return dict_left_primers_attribute,dha_temp,seq_uha_left_point,type3_left_primer_add_seq

#Fill in the "successful" and "dha" elements of dict_primers_whole based on the downstream homologous sequence primer design results, and obtain the first round of verification primer design template and seq_dha_right_point
def rightPrimerAttribute(ID,dict_right_primers,target_seq,dict_input_seq,screeningmarker_seq,config,plasmidseq):
    dict_right_primers_attribute={}
    dict_right_primers_attribute['ID']=ID
    dha_left_primer=dict_right_primers['PRIMER_LEFT_0_SEQUENCE']
    dha_right_primer=dict_right_primers['PRIMER_RIGHT_0_SEQUENCE']
    max_verify_1_up_ponit=int(config['max_verify_1_up_ponit'])
    max_verify_1_down_ponit=int(config['max_verify_1_down_ponit'])
    max_verify_2_down_ponit=int(config['max_verify_2_down_ponit'])
    length_homo_seq=int(config['length_homologous_sequence'])
    dha_right_primer_rev=revComp(dha_right_primer)
    dha_temp=dict_input_seq['seq_dha_max_whole'][:int(config['max_right_arm_seq_length'])]
    seq_dha_right_point=dha_temp.find(dha_right_primer_rev)+len(dha_right_primer)-1
    #rev_dha_right_primer_add_seq=plasmidseq[:20]
    dha_left_primer_add_seq=screeningmarker_seq[0-length_homo_seq:]
    seq_verify1_temp = ""
    # 处理上传质粒
    dha_right_primer = revComp(dha_right_primer_rev + plasmidseq[:length_homo_seq])
    # verify1 template
    if dict_input_seq["level"] == ">30&<90":
        # 处理insertion substitution的 30< alt <90
        seq_verify1_temp=dict_input_seq["seq_uha_max_whole"][-max_verify_1_up_ponit:]+dict_input_seq["seq_dr_uha"]+ screeningmarker_seq + dict_input_seq["seq_dr_dha"] + dict_input_seq["seq_dha_max_whole"][:max_verify_1_down_ponit]
        dict_right_primers_attribute['seq_verify1_temp']=seq_verify1_temp
    elif dict_input_seq["level"] == "del":
        # 处理deletion
        seq_verify1_temp=dict_input_seq["seq_uha_max_whole"][-max_verify_1_up_ponit:]+ screeningmarker_seq + dict_input_seq["seq_dr"] + dict_input_seq["seq_dha_max_whole"][:max_verify_1_down_ponit]
        dict_right_primers_attribute['seq_verify1_temp']=seq_verify1_temp
    elif dict_input_seq["level"] == ">90":
        # 处理insertion substitution的 alt >90
        seq_verify1_temp=dict_input_seq["seq_uha_max_whole"][-max_verify_1_up_ponit:] + dict_input_seq["seq_altered"] + screeningmarker_seq  + dict_input_seq["seq_dr"]+ dict_input_seq["seq_dha_max_whole"][:max_verify_1_down_ponit]
        dict_right_primers_attribute['seq_verify1_temp']=seq_verify1_temp
    elif dict_input_seq["level"] == "<30":
        # 处理insertion substitution的 alt <30
        seq_verify1_temp=dict_input_seq["seq_uha_max_whole"][-max_verify_1_up_ponit:]+ dict_input_seq["seq_altered"] + screeningmarker_seq + dict_input_seq["seq_dr"] + dict_input_seq["seq_dha_max_whole"][:max_verify_1_down_ponit]
        dict_right_primers_attribute['seq_verify1_temp']=seq_verify1_temp
    dict_right_primers_attribute['PRIMER_RIGHT_WHOLE_SEQUENCE']=dha_right_primer
    if dict_input_seq["level"] == ">30&<90":
        # 处理insertion substitution的 30< alt <90
        dha_left_primer_whole =dha_left_primer_add_seq+dict_input_seq["seq_dr_dha"]+ dha_left_primer
        dr_seq=dict_input_seq['seq_dr_dha']
    else:
        # 处理insertion substitution的 alt >90
        # 处理deletion  和 insertion substitution的 alt <= 30
        dha_left_primer_whole=dha_left_primer_add_seq+dict_input_seq["seq_dr"]+dha_left_primer
        dr_seq=dict_input_seq['seq_dr']
    dict_right_primers_attribute['PRIMER_LEFT_WHOLE_SEQUENCE']=dha_left_primer_whole
    seq_dha=dha_left_primer_add_seq+dr_seq+dha_temp[:seq_dha_right_point+1] + plasmidseq[:length_homo_seq]
    dict_right_primers_attribute['PRODUCT_SEQUENCE']=seq_dha
    dict_right_primers_attribute['PRODUCT_WHOLE_LENGTH']=len(seq_dha)
    seq_verify2_temp_dha_half=dr_seq+dha_temp[:seq_dha_right_point+max_verify_2_down_ponit+1]
    dict_right_primers_attribute['SEQ_VERIFY2_TEMP_DHA_HALF']=seq_verify2_temp_dha_half
    dict_right_primers_attribute['SEQ_DHA_WITHOIUTE_ADD']=dha_temp[:seq_dha_right_point+1]
    return dict_right_primers_attribute

def insertPrimerAttribute(ID,dict_insert_primers,dict_input_seq,screeningmarker_seq,config):
    dict_insert_primers_attribute={}
    length_homo_seq=int(config['length_homologous_sequence'])
    dict_insert_primers_attribute['ID']=ID
    insert_left_primer=dict_insert_primers['PRIMER_LEFT_0_SEQUENCE']
    insert_right_primer=dict_insert_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_insert_right_primer=revComp(insert_right_primer)
    dict_insert_primers_attribute['PRIMER_LEFT_WHOLE_SEQUENCE']=insert_left_primer
    insert_right_primer_add_seq=screeningmarker_seq[:length_homo_seq]
    insert_right_primer_whole=revComp(rev_insert_right_primer+insert_right_primer_add_seq)
    dict_insert_primers_attribute['PRIMER_RIGHT_WHOLE_SEQUENCE']=insert_right_primer_whole
    dict_insert_primers_attribute['PRIMER_LEFT_TM']=primer3.calcTm(insert_left_primer)
    dict_insert_primers_attribute['PRIMER_RIGHT_TM']=primer3.calcTm(insert_right_primer_whole)
    dict_insert_primers_attribute['PRODUCT_SEQUENCE']=dict_input_seq['seq_altered']+insert_right_primer_add_seq
    dict_insert_primers_attribute['PRIMER_PRODUCT_SIZE']=len(dict_input_seq['seq_altered']+insert_right_primer_add_seq)
    dict_insert_primers_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(insert_left_primer).tm
    dict_insert_primers_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(insert_right_primer_whole).tm
    dict_insert_primers_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(insert_left_primer,insert_right_primer_whole).tm
    return dict_insert_primers_attribute

#Fill in the "successful" "verify1" element of dict_primers_whole based on the second round of verification primer design results
def verify1PrimerAttribute(ID,dict_verify1_primers,seq_verify1_temp,assigned_tp1_seq=None):
    dict_verify1_primer_attribute={}
    dict_verify1_primer_attribute['ID']=ID
    verify1_left_primer=dict_verify1_primers['PRIMER_LEFT_0_SEQUENCE']
    verify1_right_primer=dict_verify1_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_verify1_right_primer=revComp(verify1_right_primer)
    seq_verify1_left_point=seq_verify1_temp.find(verify1_left_primer)
    if assigned_tp1_seq:
        rev_assigned_tp1_seq=revComp(assigned_tp1_seq)
        seq_verify1_right_point=seq_verify1_temp.find(rev_assigned_tp1_seq)+len(assigned_tp1_seq)-1
        seqverify1=seq_verify1_temp[seq_verify1_left_point:seq_verify1_right_point+1]
        dict_verify1_primer_attribute['PRIMER_RIGHT_0_SEQUENCE']=assigned_tp1_seq
        dict_verify1_primer_attribute['PRIMER_RIGHT_0_TM']=primer3.calcTm(assigned_tp1_seq)
        dict_verify1_primer_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(assigned_tp1_seq).tm
        dict_verify1_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(verify1_left_primer,assigned_tp1_seq).tm
    else:
        seq_verify1_right_point=seq_verify1_temp.find(rev_verify1_right_primer)+len(verify1_right_primer)-1
        seqverify1=seq_verify1_temp[seq_verify1_left_point:seq_verify1_right_point+1]
        dict_verify1_primer_attribute['PRIMER_RIGHT_0_SEQUENCE']=verify1_right_primer
        dict_verify1_primer_attribute['PRIMER_RIGHT_0_TM']=dict_verify1_primers['PRIMER_RIGHT_0_TM']
        dict_verify1_primer_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(verify1_right_primer).tm
        dict_verify1_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(verify1_left_primer,verify1_right_primer).tm
    dict_verify1_primer_attribute['PRIMER_LEFT_0_SEQUENCE']=verify1_left_primer
    dict_verify1_primer_attribute['PRIMER_LEFT_0_TM']=dict_verify1_primers['PRIMER_LEFT_0_TM']
    dict_verify1_primer_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(verify1_left_primer).tm
    dict_verify1_primer_attribute['PRODUCT_SEQUENCE']=seqverify1
    dict_verify1_primer_attribute['PRIMER_PRODUCT_SIZE']=len(seqverify1)
    return dict_verify1_primer_attribute

#Fill in the "successful" "verify1" and "verify2"element of dict_primers_whole based on the second round of verification primer design results
# screeningmarkerPrimerAttribute also used this function 
def verifyPrimerAttribute(ID,dict_verify_primers,seq_verify_temp):
    dict_verify_primer_attribute={}
    dict_verify_primer_attribute['ID']=ID
    verify_left_primer=dict_verify_primers['PRIMER_LEFT_0_SEQUENCE']
    verify_right_primer=dict_verify_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_verify_right_primer=revComp(verify_right_primer)
    seq_verify_left_point=seq_verify_temp.find(verify_left_primer)
    seq_verify_right_point=seq_verify_temp.find(rev_verify_right_primer)+len(verify_right_primer)-1
    seqverify=seq_verify_temp[seq_verify_left_point:seq_verify_right_point+1]
    dict_verify_primer_attribute['PRODUCT_SEQUENCE']=seqverify
    dict_verify_primer_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(verify_left_primer).tm
    dict_verify_primer_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(verify_right_primer).tm
    dict_verify_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(verify_left_primer,verify_right_primer).tm
    dict_verify_primer_attribute['PRIMER_LEFT_0_SEQUENCE']=dict_verify_primers['PRIMER_LEFT_0_SEQUENCE']
    dict_verify_primer_attribute['PRIMER_RIGHT_0_SEQUENCE']=dict_verify_primers['PRIMER_RIGHT_0_SEQUENCE']
    dict_verify_primer_attribute['PRIMER_LEFT_0_TM']=dict_verify_primers['PRIMER_LEFT_0_TM']
    dict_verify_primer_attribute['PRIMER_RIGHT_0_TM']=dict_verify_primers['PRIMER_RIGHT_0_TM']
    dict_verify_primer_attribute['PRIMER_PRODUCT_SIZE']=dict_verify_primers['PRIMER_PAIR_0_PRODUCT_SIZE']
    return dict_verify_primer_attribute

#Get the file containing the primers submitted to the sequence synthesis company
def primers_submitted_output(dict_uha,dict_dha,dict_screening_marker,ofile,dict_insert={},dict_verify1={},dict_verify2=None):
    # with open(output_dir,"w") as ofile:
    ofile.write(dict_uha['ID'] + '-1' + '\t'+ dict_uha['PRIMER_LEFT_WHOLE_SEQUENCE'] + '\n')
    ofile.write(dict_uha['ID'] + '-2' + '\t'+ dict_uha['PRIMER_RIGHT_WHOLE_SEQUENCE'] + '\n')
    ofile.write(dict_dha['ID'] + '-3' + '\t'+ dict_dha['PRIMER_LEFT_WHOLE_SEQUENCE'] + '\n')
    ofile.write(dict_dha['ID'] + '-4' + '\t'+ dict_dha['PRIMER_RIGHT_WHOLE_SEQUENCE'] + '\n')
    ofile.write(dict_screening_marker['ID'] + '-7' + '\t'+ dict_screening_marker['PRIMER_LEFT_0_SEQUENCE'] + '\n')
    ofile.write(dict_screening_marker['ID'] + '-8' + '\t'+ dict_screening_marker['PRIMER_RIGHT_0_SEQUENCE'] + '\n')
    if dict_insert:
        ofile.write(dict_insert['ID'] + '-5' + '\t'+ dict_insert['PRIMER_LEFT_WHOLE_SEQUENCE'] + '\n')
        ofile.write(dict_insert['ID'] + '-6' + '\t'+ dict_insert['PRIMER_RIGHT_WHOLE_SEQUENCE'] + '\n')
    if dict_verify1:
        ofile.write('test-' + dict_verify1['ID'] + '-1' + '\t'+ dict_verify1['PRIMER_LEFT_0_SEQUENCE'] + '\n')
        ofile.write('test-' + dict_verify1['ID'] + '-2' + '\t'+ dict_verify1['PRIMER_RIGHT_0_SEQUENCE'] + '\n')
    if dict_verify2:
        ofile.write('test-' + dict_verify2['ID'] + '-3' + '\t'+ dict_verify2['PRIMER_LEFT_0_SEQUENCE'] + '\n')
        ofile.write('test-' + dict_verify2['ID'] + '-4' + '\t'+ dict_verify2['PRIMER_RIGHT_0_SEQUENCE'] + '\n')
        
#Map the target sequence to the reference genome by Blast
def blast_ha(ref_genome,workdir,blast_input_file_path):
    blast_output_file_path=workdir+'/blast_output.txt'
    ref_lib=ref_genome.split('/')[-1].split('.')[0]
    seq_length=0
    with open(blast_input_file_path,'r') as ifile:
        for line in ifile:
            if not line[0]=='>':
                seq_length += len(line)-1
                break
    if seq_length > 550:
        evalue='300'
    else:
        evalue=str(int((seq_length*0.5521-7.5856)*0.8))
    os.system("makeblastdb -in "+ref_genome+" -dbtype nucl -parse_seqids -out "+ref_lib)
    os.system("blastn -query "+blast_input_file_path+" -db "+ref_lib+" -outfmt 6 -out "+blast_output_file_path+" -evalue 1e-"+evalue+" -max_target_seqs 5 -num_threads 4")
    os.system("rm %s.n*" % ref_lib)

#Evaluate the feasibility of design with the mapping of the homologous arm to the reference genome
def blast_output_evaluate(workdir,ref_genome):
    evaluate_output_file_path=workdir+'/Evaluation_result.txt'
    with open(evaluate_output_file_path,'w') as evaluate_output:
        evaluate_output.write("ID\tWarning\n")
        dict_evaluate_output={}
        with open(workdir+'/blast_output.txt','r') as evaluate_input:
            dict_result_id = {}
            for line_result in evaluate_input:
                result_id=line_result.split('\t')[0]
                dict_result_id[result_id] = dict_result_id.get(result_id,0) + 1
            list_result_id_unmap=[]
            with open(workdir+'/blast_input.txt','r') as blast_input:
                for lines in blast_input:
                    if lines[0] =='>':
                        blast_input_id=lines[1:-1]
                        if blast_input_id in dict_result_id:
                            if dict_result_id[blast_input_id]>1:
                                evaluate_output.write(blast_input_id+'\t'+'The target sequence can map to multiple positions in the reference genome. The genome editing may be mislocated.'+'\n')
                            else:
                                continue
                        else:
                            if blast_input_id in list_result_id_unmap:
                                continue
                            else:
                                evaluate_output.write(blast_input_id+'\t'+'The target sequence can not map to the reference genome. Please check them.'+'\n')
                                list_result_id_unmap.append(blast_input_id)
        for key3 in dict_evaluate_output:
            evaluate_output.write(key3+'\t'+dict_evaluate_output[key3]+'\n')
    evaluate_output_dir=evaluate_output_file_path.replace('.txt','.xlsx')
    pd.read_table(evaluate_output_file_path, index_col=0).to_excel(evaluate_output_dir)

def genOutputFile(dict_primers_whole,workdir):
    primer_succceed_dir=workdir+"/Design_results.xlsx"
    primer_failed_dir=workdir+"/Failed_task.xlsx"
    writer_succceed = pd.ExcelWriter(primer_succceed_dir)
    writer_failed = pd.ExcelWriter(primer_failed_dir)
    for kind in dict_primers_whole:
        if kind == 'successful':
            for sub_kind in dict_primers_whole[kind]:
                pd.read_json(json.dumps(dict_primers_whole[kind][sub_kind])).to_excel(writer_succceed, "Primers_for_"+sub_kind)
            writer_succceed.save()
        else:
            for sub_kind in dict_primers_whole[kind]:
                pd.read_json(json.dumps(dict_primers_whole[kind][sub_kind])).to_excel(writer_failed, "Failed_task_for_"+sub_kind)
            writer_failed.save()

def generate_visualize_file(key1,dict_input_seq,plasmidseq,screeningmarker_seq,config,dict_left_primer,dict_right_primer,dict_uha,dict_dha,dict_screeningmarker,dict_verify1,dict_verify2,fsave,workdir,dict_insert={}):
    """
    Arguments:
        key1[str]: seq id
        dict_input_seq[dict]: {
                        "seq_uha_max_whole":,
                        "seq_dha_max_whole":,
                        "seq_altered":,
                        "type":,
                        "ref":,
                        "level":,
                        "strand":,
                        "mutation_pos_index":,
                        "geneid":,
                    }
        plasmidseq[str]: plasmid sequence
        screeningmarker_seq[str]: screeningmarker_seq sequence
        config[dict]: points information
        dict_left_primer[dict]:
        dict_right_primer[dict]:
        dict_uha[dict]: uha design information
        dict_dha[dict]: dha design information
        dict_screeningmarker[dict]: screeningmarker information
        dict_verify1[dict]: verify1 design information
        dict_verify2[dict]: verify2 design information
        fsave[handler]: file handler
        workdir[str]: output dir
        dict_insert[dict]: insert primer design information
    Return: table,json
        table format:
            seqid,pos,ref,alt,primer-1,primer-2,primer-3,primer-4,test-primer-1,test-primer-2,test-primer-3,
            s3-json
        json format:     
    """
    # variables offset bp
    v_offset = 500
    
    # unpack dict_input_seq dict
    mutation_type = dict_input_seq["type"]
    uha_max = dict_input_seq["seq_uha_max_whole"]
    dha_max = dict_input_seq["seq_dha_max_whole"]
    ref = dict_input_seq["ref"]
    alt = dict_input_seq["seq_altered"]
    strand = dict_input_seq["strand"]
    geneid = dict_input_seq["geneid"]
    mutation_pos_index = dict_input_seq["mutation_pos_index"]

    dr_seq_length=int(config['dr_seq_length'])
    # visualize json
    visualize_list = []

    # 0. target seq
    targetseq = uha_max + dha_max if mutation_type == "deletion" else uha_max + alt + dha_max    
    # 1.original  seq visualize
    tmp = {}
    tmp["name"] = "target manipulation site"
    tmp["detail"] = "The position of target manipulation site (the sequence before the manipulation) on the reference genome is visualized."
    if mutation_type in ['deletion','substitution']:
        refseq = uha_max + ref + dha_max
        tmp["params"] = [
            {
                "name":"target_seq_mutation_site",
                "start":str(len(uha_max)),
                "end":str(len(uha_max)+len(ref)),
                "mutation_type":mutation_type,
                "alt":alt,
            }
        ]
    else:
        refseq = uha_max + dha_max
        tmp["params"] = [
            {
                "name":"target_seq_mutation_site",
                "start":str(len(uha_max)),
                "end":str(len(uha_max)+1),
                "mutation_type":mutation_type,
                "alt":uha_max[-1] + alt,
            }
        ]
    tmp["seq"] = refseq
    
    visualize_list.append(tmp)

    # 2.plasmid intergrated seq visualize
    # or linear sequence
    seq_uha_left_point = targetseq.find(
        dict_left_primer["PRIMER_LEFT_0_SEQUENCE"]
        )
    seq_dha_right_point = targetseq.find(
        revComp(
            dict_right_primer["PRIMER_RIGHT_0_SEQUENCE"]
            )
        ) + len(
            dict_right_primer["PRIMER_RIGHT_0_SEQUENCE"])-1
    uha = targetseq[seq_uha_left_point:len(uha_max)]
    dha = targetseq[len(uha_max)+ len(alt):seq_dha_right_point+1]

    # 根据level组装模板
    if dict_input_seq["level"] == "del":
        intergrated_seq = uha + screeningmarker_seq + dict_input_seq["seq_dr"] + dha
    elif dict_input_seq["level"] == "<30":
        intergrated_seq = uha + dict_input_seq["seq_altered"] + screeningmarker_seq + dict_input_seq["seq_dr"] + dha
    elif dict_input_seq["level"] == ">30&<90":
        intergrated_seq = uha + dict_input_seq["seq_dr_uha"]+ screeningmarker_seq + dict_input_seq["seq_dr_dha"] + dha
    elif dict_input_seq["level"] == ">90":
        intergrated_seq = uha + dict_input_seq["seq_altered"] + screeningmarker_seq + dict_input_seq["seq_dr"] + dha
    # 如果有质粒，两端添加质粒
    if plasmidseq:
        plasmid_intergrated_seq = plasmidseq[20:] + intergrated_seq + plasmidseq[:20]
    else:
        # 没有线性序列
        plasmid_intergrated_seq = intergrated_seq
    
    ## 2.1 uha
    tmp ={}
    tmp["name"] = "plasmid" if plasmidseq else "Fragment for 1st crossover"
    tmp["detail"] = "Constructed recombinant plasmid. The upstream homologous arm (UHA), the downstream homologous arm (DHA), screening marker and primer-1/2, primer-3/4, primer-5/6 (for long fragment insertion or substitution) primer-7/8 are visualized." if plasmidseq else "Assembled PCR product. The upstream homologous arm (UHA), the downstream homologous arm (DHA), screening marker and primer-1/2, primer-3/4, primer-5/6 (for long fragment insertion or substitution) primer-7/8 are visualized."
    tmp["seq"] = plasmid_intergrated_seq
    tmp["params"] = [
        {
            'name':"uha",
            'start':str(plasmid_intergrated_seq.find(uha)),
            'end':str(plasmid_intergrated_seq.find(uha) + len(uha)),
            'direction':1,
            'color':'red',
        }
    ]
    
    ## 2.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(plasmid_intergrated_seq.find(dha)),
            "end":str(plasmid_intergrated_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )

    # screeningmarker
    tmp["params"].append(
        {
            "name":"screeningmarker",
            "start":str(plasmid_intergrated_seq.find(screeningmarker_seq)),
            "end":str(plasmid_intergrated_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
            'direction':1,
            'color':'pink',
        }
    )
    # dr seq
    # 1. front dr
    tmp["params"].append(
        {
            "name":"seq_dr",
            "start":str(plasmid_intergrated_seq.find(screeningmarker_seq)-dr_seq_length),
            "end":str(plasmid_intergrated_seq.find(screeningmarker_seq)),
            'direction':1,
            'color':'brown',
        }
    )
    # 2. backend dr
    tmp["params"].append(
        {
            "name":"seq_dr",
            "start":str(plasmid_intergrated_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
            "end":str(plasmid_intergrated_seq.find(screeningmarker_seq) + len(screeningmarker_seq) + dr_seq_length),
            'direction':1,
            'color':'brown',
        }
    )
    # Inserted fragment
    if mutation_type in ['substitution','insertion']:
        tmp["params"].append(
        {
            "name":"Inserted fragment",
            "start":str(plasmid_intergrated_seq.find(uha) + len(uha)),
            "end":str(plasmid_intergrated_seq.find(screeningmarker_seq)),
            'direction':1,
            'color':'green',
        }
    )
    ## 2.3 primer1
    primer1 = dict_uha['PRIMER_LEFT_WHOLE_SEQUENCE']
    tmp["params"].append(
        {
            "name":"primer-1",
            "start":str(plasmid_intergrated_seq.find(primer1)),
            "end":str(plasmid_intergrated_seq.find(primer1) + len(primer1)),
            'direction':1,
            'color':'orange',
        }
    )

    
    ## 2.3 primer2
    primer2 = dict_uha['PRIMER_RIGHT_WHOLE_SEQUENCE']
    primer2_rev = revComp(primer2)
    tmp["params"].append(
        {
            "name":"primer-2",
            "start":str(plasmid_intergrated_seq.find(primer2_rev)),
            "end":str(plasmid_intergrated_seq.find(primer2_rev) + len(primer2_rev)),
            'direction':-1,
            'color':'orange',
        }
    )
    
    ## 2.4 primer3
    primer3 = dict_dha['PRIMER_LEFT_WHOLE_SEQUENCE']
    tmp["params"].append(
        {
            "name":"primer-3",
            "start":str(plasmid_intergrated_seq.find(primer3)),
            "end":str(plasmid_intergrated_seq.find(primer3) + len(primer3)),
            'direction':1,
            'color':'blue',
        }
    )
    
    ## 2.4 primer4
    primer4 = dict_dha['PRIMER_RIGHT_WHOLE_SEQUENCE']
    primer4_rev = revComp(primer4)
    tmp["params"].append(
        {
            "name":"primer-4",
            "start":str(plasmid_intergrated_seq.find(primer4_rev)),
            "end":str(plasmid_intergrated_seq.find(primer4_rev) + len(primer4_rev)),
            'direction':-1,
            'color':'blue',
        }
    )
    
    # 2.5 primer5
    if dict_insert:
        primer5 = dict_insert['PRIMER_LEFT_WHOLE_SEQUENCE']
        tmp["params"].append(
            {
                "name":"primer-5",
                # "start":str(plasmid_intergrated_seq.find(primer5)),
                "start":str(plasmid_intergrated_seq.find(uha) + len(uha)),
                # "end":str(plasmid_intergrated_seq.find(primer5) + len(primer5)),
                "end":str(plasmid_intergrated_seq.find(uha) + len(uha) + len(primer5)),
                'direction':1,
                'color':'#0dcaf0',
            }
        )

        # 2.6 primer6
        primer6 = dict_insert['PRIMER_RIGHT_WHOLE_SEQUENCE']
        primer6_rev = revComp(primer6)
        tmp["params"].append(
            {
                "name":"primer-6",
                "start":str(plasmid_intergrated_seq.find(primer6_rev)),
                "end":str(plasmid_intergrated_seq.find(primer6_rev) + len(primer6_rev)),
                'direction':-1,
                'color':'#0dcaf0',
            }
        )
    else:
        primer5 = ""
        primer6 = ""

    # p7 p8 for screeningmarker
    # 2.7
    primer7 = dict_screeningmarker['PRIMER_LEFT_0_SEQUENCE']
    tmp["params"].append(
        {
            "name":"primer-7",
            "start":str(plasmid_intergrated_seq.find(primer7)),
            "end":str(plasmid_intergrated_seq.find(primer7) + len(primer7)),
            'direction':1,
            'color':'purple',
        }
    )
    # 2.8
    primer8 = dict_screeningmarker['PRIMER_RIGHT_0_SEQUENCE']
    primer8_rev = revComp(primer8)
    tmp["params"].append(
        {
            "name":"primer-8",
            "start":str(plasmid_intergrated_seq.find(primer8_rev)),
            "end":str(plasmid_intergrated_seq.find(primer8_rev) + len(primer8_rev)),
            'direction':-1,
            'color':'purple',
        }
    )
    visualize_list.append(tmp)

    # 3. verify1 visualize
    # 有质粒，才有v1
    if plasmidseq:
        v1_seq = dict_input_seq["uha_upstream"]+ targetseq[:seq_uha_left_point] + intergrated_seq + plasmidseq[
            :int(config["max_verify_1_down_ponit"]) + v_offset
        ]
        
        ## 3.1 uha
        tmp = {}
        tmp["name"]="Verify1"
        tmp["detail"]="Verification of the 1st-round of crossover and isolation. The upstream homologous arm (UHA), the downstream homologous arm (DHA) screening marker and test-primer-1/2 are visualized."
        tmp["seq"] = v1_seq
        tmp["params"] = [
            {
                "name":"uha",
                "start":str(v1_seq.find(uha)),
                "end":str(v1_seq.find(uha) + len(uha)),
                'direction':1,
                'color':'red',
            }
        ]
        
        ## 3.2 dha
        tmp["params"].append(
            {
                "name":"dha",
                "start":str(v1_seq.find(dha)),
                "end":str(v1_seq.find(dha) + len(dha)),
                'direction':1,
                'color':'red',
            }
        )

        # screeningmarker
        tmp["params"].append(
            {
                "name":"screeningmarker",
                "start":str(v1_seq.find(screeningmarker_seq)),
                "end":str(v1_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
                'direction':1,
                'color':'pink',
            }
        )
        # dr seq
        # 1. front dr
        tmp["params"].append(
            {
                "name":"seq_dr",
                "start":str(v1_seq.find(screeningmarker_seq)-dr_seq_length),
                "end":str(v1_seq.find(screeningmarker_seq)),
                'direction':1,
                'color':'brown',
            }
        )
        # 2. backend dr
        tmp["params"].append(
            {
                "name":"seq_dr",
                "start":str(v1_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
                "end":str(v1_seq.find(screeningmarker_seq) + len(screeningmarker_seq) + dr_seq_length),
                'direction':1,
                'color':'brown',
            }
        )
        ## Inserted fragment
        if mutation_type in ['substitution','insertion']:
            tmp["params"].append(
            {
                "name":"Inserted fragment",
                "start":str(v1_seq.find(uha) + len(uha)),
                "end":str(v1_seq.find(screeningmarker_seq)),
                'direction':1,
                'color':'green',
            }
        )
        ## 3.3 v1 primer1
        v1_primer1 = dict_verify1['PRIMER_LEFT_0_SEQUENCE']
        tmp["params"].append(
            {
                "name":"test-primer-1",
                "start":str(v1_seq.find(v1_primer1)),
                "end":str(v1_seq.find(v1_primer1) + len(v1_primer1)),
                'direction':1,
                'color':'orange',
            }
        )
        
        
        ## 3.4 v1 primer2
        v1_primer2 = dict_verify1['PRIMER_RIGHT_0_SEQUENCE']
        v1_primer2_rev =revComp(v1_primer2)
        tmp["params"].append(
            {
                "name":"test-primer-2",
                "start":str(v1_seq.find(v1_primer2_rev)),
                "end":str(v1_seq.find(v1_primer2_rev) + len(v1_primer2_rev)),
                'direction':-1,
                'color':'orange',
            }
        )
        visualize_list.append(tmp)
    
    # 4. verify2 visualize
    # v2_seq = dict_input_seq["uha_upstream"]+ targetseq[:seq_uha_left_point] + uha + alt + dha + targetseq[seq_dha_right_point + 1 :]  + dict_input_seq["dha_downstream"]
    v2_seq = dict_input_seq["uha_upstream"]+ targetseq + dict_input_seq["dha_downstream"]
    
    ## 4.1 uha
    tmp ={}
    tmp["name"] = "Verify2"
    tmp["detail"] = "Verification of the 2nd-round of crossover and isolation. The upstream homologous arm (UHA), the downstream homologous arm (DHA) and test-primer-3/4 are visualized."
    tmp["seq"] = v2_seq
    tmp["params"] = [
        {
            "name":"uha",
            "start":str(v2_seq.find(uha)),
            "end":str(v2_seq.find(uha) + len(uha)),
            'direction':1,
            'color':'red',
        }
    ]
    
    ## 4.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(v2_seq.find(dha)),
            "end":str(v2_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )
    ## Inserted fragment
    if mutation_type in ['substitution','insertion']:
        tmp["params"].append(
        {
            "name":"Inserted fragment",
            "start":str(v2_seq.find(uha) + len(uha)),
            "end":str(v2_seq.find(dha)),
            'direction':1,
            'color':'green',
        }
    )
    ## 4.3 v2 primer1
    v2_primer1 = dict_verify2['PRIMER_LEFT_0_SEQUENCE']
    tmp["params"].append(
        {
            "name":"test-primer-3",
            "start":str(v2_seq.find(v2_primer1)),
            "end":str(v2_seq.find(v2_primer1) + len(v2_primer1)),
            'direction':1,
            'color':'orange',
        }
    )
    
    ## 4.4 v2 primer2
    v2_primer2 = dict_verify2['PRIMER_RIGHT_0_SEQUENCE']
    v2_primer2_rev =revComp(v2_primer2)
    tmp["params"].append(
        {
            "name":"test-primer-4",
            "start":str(v2_seq.find(v2_primer2_rev)),
            "end":str(v2_seq.find(v2_primer2_rev) + len(v2_primer2_rev)),
            'direction':-1,
            'color':'orange',
        }
    )
    visualize_list.append(tmp)
    
    # write to json
    json_file = os.path.join(
        workdir,
        '%s.json' % key1
    )
    with open(json_file,"w") as f:
        json.dump(visualize_list,f,indent=4)
    
    # write to table file , fsave handler
    linelist = [
        key1,
        ref,
        alt,
        mutation_type,
        dict_uha['PRIMER_LEFT_WHOLE_SEQUENCE'],
        dict_uha['PRIMER_RIGHT_WHOLE_SEQUENCE'],
        dict_dha['PRIMER_LEFT_WHOLE_SEQUENCE'],
        dict_dha['PRIMER_RIGHT_WHOLE_SEQUENCE'],
        primer5,
        primer6,
        primer7,
        primer8,
        dict_verify1['PRIMER_LEFT_0_SEQUENCE'] if plasmidseq else "",
        dict_verify1['PRIMER_RIGHT_0_SEQUENCE'] if plasmidseq else "",
        dict_verify2['PRIMER_LEFT_0_SEQUENCE'],
        dict_verify2['PRIMER_RIGHT_0_SEQUENCE'],
        strand,
        geneid,
        str(mutation_pos_index),
    ]
    fsave.write(
        '\t'.join(linelist) + "\n"
    )

#Design process
def design_process(input_file_path,screeningmarker_file_path,workdir,ref_genome,config,plasmid_file_path,assigned_tp1_seq=None,screening_maker_removal=None):
    dict_input_seq=input_to_primer_template(input_file_path,ref_genome,config,workdir)
    if isinstance(dict_input_seq,str):
        #print('*********** error ***************')
        #print(dict_input_seq)
        return dict_input_seq
    dict_primers_whole={'successful':{'uha':[],'dha':[],'verify1':[],'verify2':[],'inserted':[],'screening_marker':[]},'failed':{'uha':[],'dha':[],'verify1':[],'verify2':[],'inserted':[],'screening_marker':[]}}
    screeningmarker_seq=linearSeq(screeningmarker_file_path)
    plasmidseq=linearSeq(plasmid_file_path)
    failed_sheet_title=['PRIMER_LEFT_EXPLAIN','PRIMER_RIGHT_EXPLAIN','PRIMER_PAIR_EXPLAIN','PRIMER_LEFT_NUM_RETURNED','PRIMER_RIGHT_NUM_RETURNED','PRIMER_INTERNAL_NUM_RETURNED','PRIMER_PAIR_NUM_RETURNED']
    # 读写文件
    blast_input_file_path=workdir+'/blast_input.txt'
    primer_order_file_path=workdir+'/Primer_order.txt'
    visualize_file_path=workdir+'/Visualize_order.tsv'
    blast_input_file =  open(blast_input_file_path,'w')
    primer_order_fsave =  open(primer_order_file_path,'w')
    visualize_fsave =  open(visualize_file_path,'w')
    # write visualize file header
    headers = ['seqid','ref','alt',"manipulation","primer-1","primer-2","primer-3","primer-4","primer-5","primer-6","primer-7","primer-8","test-primer-1","test-primer-2","test-primer-3","test-primer-4","strand","geneid","mutation_pos_index"]
    visualize_fsave.write(
        '\t'.join(headers) + "\n"
    )
    for key1 in dict_input_seq:
        dict_inserted_primers_attribute ={}
        dictscreeningmarkerprimers=primerDesign(key1,screeningmarker_seq,config,"screening_marker")
        if len(dictscreeningmarkerprimers)<10:
            dict_screening_marker_failed_output={k:dictscreeningmarkerprimers[k] for k in dictscreeningmarkerprimers.keys() if k in failed_sheet_title}
            dict_screening_marker_failed_output['ID']=key1
            dict_primers_whole['failed']['screening_marker'].append(dict_screening_marker_failed_output)
        else:
            dict_screening_marker_primers_attribute = verifyPrimerAttribute(key1,dictscreeningmarkerprimers,screeningmarker_seq)
            dict_primers_whole['successful']['screening_marker'].append(dict_screening_marker_primers_attribute)
            target_seq=dict_input_seq[key1]['seq_uha_max_whole']+dict_input_seq[key1]['seq_altered']+dict_input_seq[key1]['seq_dha_max_whole']
            uha_temp=dict_input_seq[key1]['seq_uha_max_whole'][len(dict_input_seq[key1]['seq_uha_max_whole'])-int(config['max_left_arm_seq_length']):]
            dictleftprimers=primerDesign(key1,uha_temp,config,"left_arm")
            if len(dictleftprimers)<10:
                dict_uha_failed_output={k:dictleftprimers[k] for k in dictleftprimers.keys() if k in failed_sheet_title}
                dict_uha_failed_output['ID']=key1
                dict_primers_whole['failed']['uha'].append(dict_uha_failed_output)
                print("error")
            else:
                dict_left_primers_attribute = leftPrimerAttribute(key1,dictleftprimers,target_seq,dict_input_seq[key1],screeningmarker_seq,config,plasmidseq)
                dict_primers_whole['successful']['uha'].append(dict_left_primers_attribute)
                dha_temp=dict_input_seq[key1]['seq_dha_max_whole'][:int(config['max_right_arm_seq_length'])]
                dictrightprimers=primerDesign(key1,dha_temp,config,"right_arm")
                if len(dictrightprimers)<10:
                    dict_dha_failed_output={k:dictrightprimers[k] for k in dictrightprimers.keys() if k in failed_sheet_title}
                    dict_dha_failed_output['ID']=key1
                    dict_primers_whole['failed']['dha'].append(dict_dha_failed_output)
                    print("error")
                else:
                    dict_right_primers_attribute = rightPrimerAttribute(key1,dictrightprimers,target_seq,dict_input_seq[key1],screeningmarker_seq,config,plasmidseq)
                    dict_primers_whole['successful']['dha'].append(dict_right_primers_attribute)
                    #seq_verify1_temp=rightPrimerAttribute(key1,dictrightprimers,target_seq,dict_input_seq[key1],plasmidseq,config,seq_uha_left_point,"type3")[1]
                    #seq_dha_right_point=rightPrimerAttribute(key1,dictrightprimers,target_seq,dict_input_seq[key1],screeningmarker_seq,config)[1]
                    homo_arm_seq_1=dict_left_primers_attribute['SEQ_UHA_WITHOIUTE_ADD']
                    homo_arm_seq_2=dict_right_primers_attribute['SEQ_DHA_WITHOIUTE_ADD']
                    blast_input_file.write('>'+key1+'_UHA'+'\n'+homo_arm_seq_1+'\n')
                    blast_input_file.write('>'+key1+'_DHA'+'\n'+homo_arm_seq_2+'\n')
                    #dictverify1primers=primerDesign(key1,seq_verify1_temp,config,"verify_1","type3")
                    seq_verify1_temp = dict_right_primers_attribute['seq_verify1_temp']
                    dictverify1primers=primerDesign(key1,seq_verify1_temp,config,"verify_1")
                    if len(dictverify1primers)<10:
                        dict_verify1_failed_output={k:dictverify1primers[k] for k in dictverify1primers.keys() if k in failed_sheet_title}
                        dict_verify1_failed_output['ID']=key1
                        dict_primers_whole['failed']['verify1'].append(dict_verify1_failed_output)
                        print("error")
                    else:
                        if assigned_tp1_seq:
                            dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp,assigned_tp1_seq)
                        else:
                            dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp)
                        dict_primers_whole['successful']['verify1'].append(dict_verify1_primers_attribute)
                    if dict_input_seq[key1]["level"] == ">90":
                        inserted_temp = dict_input_seq[key1]['seq_altered']
                        dictinsertedprimers=primerDesign(key1,inserted_temp,config,"insert")
                        if len(dictinsertedprimers)<10:
                            dict_inserted_failed_output={k:dictinsertedprimers[k] for k in dictinsertedprimers.keys() if k in failed_sheet_title}
                            dict_inserted_failed_output['ID']=key1
                            dict_primers_whole['failed']['inserted'].append(dict_inserted_failed_output)
                            print("error")
                        else:
                            dict_inserted_primers_attribute = insertPrimerAttribute(key1,dictinsertedprimers,dict_input_seq[key1],screeningmarker_seq,config)
                            dict_primers_whole['successful']['inserted'].append(dict_inserted_primers_attribute)
                            seq_verify2_temp = dict_left_primers_attribute['SEQ_VERIFY2_TEMP_UHA_HALF'] + dict_right_primers_attribute['SEQ_VERIFY2_TEMP_DHA_HALF']
                            dictverify2primers=primerDesign(key1,seq_verify2_temp,config,"verify_2")
                            if len(dictverify2primers)<10:
                                dict_verify2_failed_output={k:dictverify2primers[k] for k in dictverify2primers.keys() if k in failed_sheet_title}
                                dict_verify2_failed_output['ID']=key1
                                dict_primers_whole['failed']['verify2'].append(dict_verify2_failed_output)
                                print("error")
                            else:
                                dict_verify2_primers_attribute = verifyPrimerAttribute(key1,dictverify2primers,seq_verify2_temp)
                                dict_primers_whole['successful']['verify2'].append(dict_verify2_primers_attribute)
                                primers_submitted_output(dict_left_primers_attribute,dict_right_primers_attribute,dict_screening_marker_primers_attribute,primer_order_fsave,dict_inserted_primers_attribute,dict_verify1_primers_attribute,dict_verify2_primers_attribute)
                    else:
                        if dict_input_seq[key1]['type']=="deletion" and screening_maker_removal=="No":
                            primers_submitted_output(dict_left_primers_attribute,dict_right_primers_attribute,dict_screening_marker_primers_attribute,primer_order_fsave,dict_verify1=dict_verify1_primers_attribute)
                        else:
                            seq_verify2_temp = dict_left_primers_attribute['SEQ_VERIFY2_TEMP_UHA_HALF'] + dict_right_primers_attribute['SEQ_VERIFY2_TEMP_DHA_HALF']
                            dictverify2primers=primerDesign(key1,seq_verify2_temp,config,"verify_2")
                            if len(dictverify2primers)<10:
                                dict_verify2_failed_output={k:dictverify2primers[k] for k in dictverify2primers.keys() if k in failed_sheet_title}
                                dict_verify2_failed_output['ID']=key1
                                dict_primers_whole['failed']['verify2'].append(dict_verify2_failed_output)
                                print("error")
                            else:
                                dict_verify2_primers_attribute = verifyPrimerAttribute(key1,dictverify2primers,seq_verify2_temp)
                                dict_primers_whole['successful']['verify2'].append(dict_verify2_primers_attribute)
                            #primers_submitted_output(leftPrimerAttribute(key1,dictleftprimers,target_seq,dict_input_seq[key1],screeningmarker_seq,config)[0],rightPrimerAttribute(key1,dictrightprimers,target_seq,dict_input_seq[key1],screeningmarker_seq,seq_uha_left_point,config)[0],verify2PrimerAttribute(key1,dictverify2primers,seq_verify2_temp),primer_order_fsave)
                                primers_submitted_output(dict_left_primers_attribute,dict_right_primers_attribute,dict_screening_marker_primers_attribute,primer_order_fsave,dict_verify1=dict_verify1_primers_attribute,dict_verify2=dict_verify2_primers_attribute)
                    # visualize
                    generate_visualize_file(
                        key1,
                        dict_input_seq[key1],
                        plasmidseq,
                        screeningmarker_seq,
                        config,
                        dictleftprimers,
                        dictrightprimers,
                        dict_left_primers_attribute,
                        dict_right_primers_attribute,
                        dict_screening_marker_primers_attribute,
                        dict_verify1_primers_attribute,
                        dict_verify2_primers_attribute,
                        visualize_fsave,
                        workdir,
                        dict_inserted_primers_attribute
                        )
    #print(dict_primers_whole['successful']['dha'])
    # save result files
    blast_input_file.close()
    primer_order_fsave.close()
    visualize_fsave.close()
    primer_order_dir=primer_order_file_path.replace('.txt','.xlsx')
    pd.read_table(primer_order_file_path, index_col=0).to_excel(primer_order_dir)
    blast_ha(ref_genome,workdir,blast_input_file_path)
    blast_output_evaluate(workdir,ref_genome)
    genOutputFile(dict_primers_whole,workdir)
    return 'success'

if __name__ == "__main__":
    import time
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', help='input sequences file', required=True)
    parser.add_argument('--plasmid', '-p', help='plasmid file', required=True)
    parser.add_argument('--screeningmarker', '-s', help='screening marker seq file', required=True)
    parser.add_argument('--ref', '-r', help='reference genome', required=True)
    parser.add_argument('--dir', '-d', help='outputdir', required=True)
    parser.add_argument('--conf', '-c', help='config', required=True)
    parser.add_argument('--ustp', '-u', help='user specific test primer', required=False)
    parser.add_argument('--smr', '-m', help='screening maker removal', required=False)
    arguments = parser.parse_args()
    input_file_path = arguments.input
    plasmid_file_path = arguments.plasmid
    screeningmarker_file_path = arguments.screeningmarker
    ref_genome = arguments.ref
    workdir = arguments.dir
    point_dict = conf_read(arguments.conf)
    user_specific_test_primer = arguments.ustp
    screening_maker_removal = arguments.smr
    time1=time.time()
    design_process(input_file_path,screeningmarker_file_path,workdir,ref_genome,point_dict,plasmid_file_path,user_specific_test_primer,screening_maker_removal)
    time2=time.time()
    print(time2-time1)