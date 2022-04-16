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
                    "reverse_mapped":1 if (int(start) > int(end) and int(float(identity)) == 100 and allength==fasta_length_dict[key]) else 0,
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
                "seq_dha_whole":"",
                "seq_altered":"",
                "seq_dr":"",
                "type":"",   # [substitution,deletion,insertion]
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_whole  down 100bp sequence,
            },
            "key2":{
                "seq_uha_max_whole":"",
                "seq_dha_whole":"",
                "seq_altered":"",
                "seq_dr_uha":"",
                "seq_dr_dha":"",
                "type":"",
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_whole  down 100bp sequence,
            }
        }
    """
    left_arm_seq_length=int(config['left_arm_seq_length'])
    right_arm_seq_length=int(config['right_arm_seq_length'])
    # 如果未上传质粒，max_verify_1_up_ponit max_verify_2_down_ponit为0
    max_verify_up_ponit=int(config['max_verify_2_up_ponit'])
    max_verify_down_ponit=int(config['max_verify_2_down_ponit'])
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
                elif blast_search_dict[mun_id]["reverse_mapped"]:
                    # 反义链
                    error_message = "The upstream sequence of " + mun_id + " can be mapped to the antisense strand, %s. Please rightly prepare input file for target manipulation as the example of 2,3-BD." % blast_search_dict[mun_id]["description"]
                    print(error_message)
                    return error_message
                elif blast_search_dict[mun_id]["unique_mapped"] == 1:
                    # 开始突变的碱基在genome上的索引
                    record = str(record_dict[blast_search_dict[mun_id]["chrom"]].seq)
                    upstream_start_index = int(blast_search_dict[mun_id]["start"])-1
                    mutation_pos_index = upstream_start_index + len(upstream)
                    mutation_pos_ref = record[mutation_pos_index].upper()
                    # length 
                    if mutation_pos_index - left_arm_seq_length - max_verify_up_ponit < 0:
                        error_message = "The length of upstream sequence of manipulation site of " + mun_id + " must be larger than sum of 'Max Length of UHA' and 'Max Length of UIS'."
                        print(error_message)
                        return error_message
                    else:
                        if mutation_type not in ["insertion","deletion","substitution"]:
                            error_message = "The target manipulation type of " + mun_id + " must be equal to 'insertion,substitution or deletion', Please rightly prepare input file for target manipulation as the example of 2,3-BD."
                            print(error_message)
                            return  error_message
                        elif mutation_type == "insertion":
                            ref = ""
                        elif mutation_type == "deletion":
                            alt = ""
                        genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]
                        if genome_ref.upper() == ref.upper():
                            primer_template[mun_id] = {
                                "seq_uha":str(record[
                                    mutation_pos_index - left_arm_seq_length:mutation_pos_index
                                    ]),
                                "seq_uha_verify":str(record[
                                    mutation_pos_index - left_arm_seq_length - max_verify_up_ponit:mutation_pos_index
                                    ]),
                                "seq_dha":str(record[
                                    mutation_pos_index + len(ref)
                                    : mutation_pos_index + len(ref) + right_arm_seq_length
                                    ]),
                                "seq_dha_verify":str(record[
                                    mutation_pos_index + len(ref)
                                    : mutation_pos_index + len(ref) + right_arm_seq_length + max_verify_down_ponit
                                    ]),
                                "seq_altered":alt,
                                "type":mutation_type,
                                "ref": ref,
                            }
                            if len(alt)<=length_threshold_msddsc:
                                primer_template[mun_id]["level"]="<=length_homologous_sequence"
                            else:
                                primer_template[mun_id]["level"]=">length_homologous_sequence"
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
    if stype =="verify1":
        print("************* %s process ************" %stype)
        left_length=100
        right_length=int(config["max_verify_2_down_ponit"])-int(config["min_verify_2_down_ponit"])
        region_list = [[0,left_length,seqlength-right_length-1,right_length]]
        size_range = [seqlength-left_length-right_length+36,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    elif stype =="verify2":
        print("************* %s process ************" %stype)
        left_length=int(config["max_verify_2_up_ponit"])-int(config["min_verify_2_up_ponit"])
        right_length=int(config["max_verify_2_down_ponit"])-int(config["min_verify_2_down_ponit"])
        region_list = [[0,left_length,seqlength-right_length-1,right_length]]
        size_range = [seqlength-left_length-right_length+36,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = region_list
    elif stype == "insert" or "screening_marker":
        print("************* %s process ************" %stype)
        size_range = [seqlength,seqlength]
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = [[0,25,seqlength-26,25]]
        seq_args["SEQUENCE_FORCE_LEFT_START"] = 0
        seq_args["SEQUENCE_FORCE_RIGHT_START"] = seqlength-1
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

def insertPrimerAttribute(ID,dict_insert_primers,dict_input_seq,input_seq):
    dict_insert_primers_attribute={}
    dict_insert_primers_attribute['ID']=ID
    insert_left_primer=dict_insert_primers['PRIMER_LEFT_0_SEQUENCE']
    insert_right_primer=dict_insert_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_insert_right_primer=revComp(insert_right_primer)
    dict_insert_primers_attribute['PRIMER_LEFT_WHOLE_SEQUENCE']=dict_input_seq['seq_uha']+insert_left_primer
    dict_insert_primers_attribute['PRIMER_RIGHT_WHOLE_SEQUENCE']=revComp(rev_insert_right_primer+dict_input_seq['seq_dha'])
    #dict_insert_primers_attribute['PRIMER_LEFT_TM']=primer3.calcTm(insert_left_primer)
    #dict_insert_primers_attribute['PRIMER_RIGHT_TM']=primer3.calcTm(insert_right_primer)
    dict_insert_primers_attribute['PRODUCT_SEQUENCE']=dict_input_seq['seq_uha']+input_seq+dict_input_seq['seq_dha']
    dict_insert_primers_attribute['PRIMER_PRODUCT_SIZE']=len(dict_input_seq['seq_uha']+input_seq+dict_input_seq['seq_dha'])
    #dict_insert_primers_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(insert_left_primer).tm
    #dict_insert_primers_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(insert_right_primer).tm
    #dict_insert_primers_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(insert_left_primer,insert_right_primer).tm
    return dict_insert_primers_attribute

def screeningmarkerPrimerAttribute(ID,dict_insert_primers,dict_input_seq,input_seq):
    dict_insert_primers_attribute={}
    dict_insert_primers_attribute['ID']=ID
    insert_left_primer=dict_insert_primers['PRIMER_LEFT_0_SEQUENCE']
    insert_right_primer=dict_insert_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_insert_right_primer=revComp(insert_right_primer)
    dict_insert_primers_attribute['PRIMER_LEFT_WHOLE_SEQUENCE']=dict_input_seq['seq_uha']+insert_left_primer
    dict_insert_primers_attribute['PRIMER_RIGHT_WHOLE_SEQUENCE']=revComp(rev_insert_right_primer+dict_input_seq['seq_dha'])
    #dict_insert_primers_attribute['PRIMER_LEFT_TM']=primer3.calcTm(insert_left_primer)
    #dict_insert_primers_attribute['PRIMER_RIGHT_TM']=primer3.calcTm(insert_right_primer)
    dict_insert_primers_attribute['PRODUCT_SEQUENCE']=dict_input_seq['seq_uha']+input_seq+dict_input_seq['seq_dha']
    dict_insert_primers_attribute['PRIMER_PRODUCT_SIZE']=len(dict_input_seq['seq_uha']+input_seq+dict_input_seq['seq_dha'])
    if dict_input_seq["level"]=="<=length_homologous_sequence":
        dict_insert_primers_attribute['OLIGONUCLEOTIDE_SEQUENCE']=dict_input_seq['seq_uha']+dict_input_seq['seq_altered']+dict_input_seq['seq_dha']
    #dict_insert_primers_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(insert_left_primer).tm
    #dict_insert_primers_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(insert_right_primer).tm
    #dict_insert_primers_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(insert_left_primer,insert_right_primer).tm
    return dict_insert_primers_attribute

#Fill in the "successful" "verify1" element of dict_primers_whole based on the second round of verification primer design results
def verify1PrimerAttribute(ID,dict_verify1_primers,seq_verify1_temp,seq_verify1_temp_add_screeningmarker=None,assigned_tp1_seq=None):
    dict_verify1_primer_attribute={}
    dict_verify1_primer_attribute['ID']=ID
    verify1_left_primer=dict_verify1_primers['PRIMER_LEFT_0_SEQUENCE']
    verify1_right_primer=dict_verify1_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_verify1_right_primer=revComp(verify1_right_primer)
    if assigned_tp1_seq:
        seq_verify1_left_point=seq_verify1_temp_add_screeningmarker.find(assigned_tp1_seq)
        print(seq_verify1_left_point)
        seq_verify1_right_point=seq_verify1_temp_add_screeningmarker.find(rev_verify1_right_primer)+len(verify1_right_primer)-1
        seqverify1=seq_verify1_temp_add_screeningmarker[seq_verify1_left_point:seq_verify1_right_point+1]
        dict_verify1_primer_attribute['PRIMER_LEFT_0_SEQUENCE']=assigned_tp1_seq
        dict_verify1_primer_attribute['PRIMER_LEFT_0_TM']=primer3.calcTm(assigned_tp1_seq)
        dict_verify1_primer_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(assigned_tp1_seq).tm
        dict_verify1_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(assigned_tp1_seq,verify1_right_primer).tm
    else:
        seq_verify1_left_point=seq_verify1_temp.find(verify1_left_primer)
        seq_verify1_right_point=seq_verify1_temp.find(rev_verify1_right_primer)+len(verify1_right_primer)-1
        seqverify1=seq_verify1_temp[seq_verify1_left_point:seq_verify1_right_point+1]
        dict_verify1_primer_attribute['PRIMER_LEFT_0_SEQUENCE']=verify1_left_primer
        dict_verify1_primer_attribute['PRIMER_LEFT_0_TM']=dict_verify1_primers['PRIMER_LEFT_0_TM']
        dict_verify1_primer_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(verify1_left_primer).tm
        dict_verify1_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(verify1_left_primer,verify1_right_primer).tm
    dict_verify1_primer_attribute['PRODUCT_SEQUENCE']=seqverify1
    dict_verify1_primer_attribute['PRIMER_RIGHT_0_SEQUENCE']=verify1_right_primer
    dict_verify1_primer_attribute['PRIMER_RIGHT_0_TM']=dict_verify1_primers['PRIMER_RIGHT_0_TM']
    dict_verify1_primer_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(verify1_right_primer).tm
    dict_verify1_primer_attribute['PRIMER_PRODUCT_SIZE']=len(seqverify1)
    return dict_verify1_primer_attribute

#Fill in the "successful" "verify2" element of dict_primers_whole based on the second round of verification primer design results
def verify2PrimerAttribute(ID,dict_verify2_primers,seq_verify2_temp):
    dict_verify2_primer_attribute={}
    dict_verify2_primer_attribute['ID']=ID
    verify2_left_primer=dict_verify2_primers['PRIMER_LEFT_0_SEQUENCE']
    verify2_right_primer=dict_verify2_primers['PRIMER_RIGHT_0_SEQUENCE']
    rev_verify2_right_primer=revComp(verify2_right_primer)
    seq_verify2_left_point=seq_verify2_temp.find(verify2_left_primer)
    seq_verify2_right_point=seq_verify2_temp.find(rev_verify2_right_primer)+len(verify2_right_primer)-1
    seqverify2=seq_verify2_temp[seq_verify2_left_point:seq_verify2_right_point+1]
    dict_verify2_primer_attribute['PRODUCT_SEQUENCE']=seqverify2
    dict_verify2_primer_attribute['PRIMER_LEFT_TM_HOMODIMER']=primer3.calcHomodimer(verify2_left_primer).tm
    dict_verify2_primer_attribute['PRIMER_RIGHT_TM_HOMODIMER']=primer3.calcHomodimer(verify2_right_primer).tm
    dict_verify2_primer_attribute['PRIMER_PAIR_TM_HETERODIMER']=primer3.calcHeterodimer(verify2_left_primer,verify2_right_primer).tm
    dict_verify2_primer_attribute['PRIMER_LEFT_0_SEQUENCE']=dict_verify2_primers['PRIMER_LEFT_0_SEQUENCE']
    dict_verify2_primer_attribute['PRIMER_RIGHT_0_SEQUENCE']=dict_verify2_primers['PRIMER_RIGHT_0_SEQUENCE']
    dict_verify2_primer_attribute['PRIMER_LEFT_0_TM']=dict_verify2_primers['PRIMER_LEFT_0_TM']
    dict_verify2_primer_attribute['PRIMER_RIGHT_0_TM']=dict_verify2_primers['PRIMER_RIGHT_0_TM']
    dict_verify2_primer_attribute['PRIMER_PRODUCT_SIZE']=dict_verify2_primers['PRIMER_PAIR_0_PRODUCT_SIZE']
    return dict_verify2_primer_attribute

#Get the file containing the primers submitted to the sequence synthesis company
def primers_submitted_output(dict_screening_marker,dict_verify1,ofile,dict_verify2=None,dict_insert={}):
    # with open(output_dir,"w") as ofile:
    ofile.write(dict_screening_marker['ID'] + '-1' + '\t'+ dict_screening_marker['PRIMER_LEFT_WHOLE_SEQUENCE'] + '\n')
    ofile.write(dict_screening_marker['ID'] + '-2' + '\t'+ dict_screening_marker['PRIMER_RIGHT_WHOLE_SEQUENCE'] + '\n')
    if dict_insert:
        ofile.write(dict_insert['ID'] + '-3' + '\t'+ dict_insert['PRIMER_LEFT_WHOLE_SEQUENCE'] + '\n')
        ofile.write(dict_insert['ID'] + '-4' + '\t'+ dict_insert['PRIMER_RIGHT_WHOLE_SEQUENCE'] + '\n')
    ofile.write('test-' + dict_verify1['ID'] + '-1' + '\t'+ dict_verify1['PRIMER_LEFT_0_SEQUENCE'] + '\n')
    ofile.write('test-' + dict_verify1['ID'] + '-2' + '\t'+ dict_verify1['PRIMER_RIGHT_0_SEQUENCE'] + '\n')
        
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

def generate_visualize_file(key1,dict_input_seq,screeningmarker_seq,config,dict_screeningmarker,dict_verify1,dict_verify2,fsave,workdir,dict_insert={}):
    """
    Arguments:
        key1[str]: seq id
        dict_input_seq[dict]: {
                        "seq_uha_max_whole":,
                        "seq_dha_whole":,
                        "seq_altered":,
                        "type":,
                        "ref":,
                    }
        screeningmarker_seq[str]: screeningmarker_seq sequence
        config[dict]: points information
        dict_screeningmarker[dict]: screeningmarker design information
        dict_verify1[dict]: verify1 design information
        dict_verify2[dict]: verify2 design information
        fsave[handler]: file handler
        workdir[str]: output dir
        dict_insert[dict]: insert primer design information
    Return: table,json
        table format:
            seqid,pos,ref,alt,primer1,primer2,primer3,primer4,test-primer-1,test-primer-2,test-primer-3,
            s3-json
        json format:     
    """
    
    #print('dict_input_seq', dict_input_seq)
    #print(key1,dict_input_seq["level"],dict_insert)
    # unpack dict_input_seq dict
    mutation_type = dict_input_seq["type"]
    uha = dict_input_seq["seq_uha"]
    dha = dict_input_seq["seq_dha"]
    ref = dict_input_seq["ref"]
    alt = dict_input_seq["seq_altered"]
    # visualize json
    visualize_list = []

    # 0. target seq
    targetseq = uha + dha if mutation_type == "deletion" else uha + alt + dha    
    # 1.original  seq visualize
    tmp = {}
    tmp["name"] = "target manipulation site"
    tmp["detail"] = "The position of target manipulation site (the sequence before the manipulation) on the reference genome is visualized."
    if mutation_type in ['deletion','substitution']:
        refseq = uha + ref + dha
        tmp["params"] = [
            {
                "name":"target_seq_mutation_site",
                "start":str(len(uha)),
                "end":str(len(uha)+len(ref)),
                "mutation_type":mutation_type,
                "alt":alt,
            }
        ]
    else:
        refseq = uha + dha
        tmp["params"] = [
            {
                "name":"target_seq_mutation_site",
                "start":str(len(uha)),
                "end":str(len(uha)+1),
                "mutation_type":mutation_type,
                "alt":uha[-1] + alt,
            }
        ]
    tmp["seq"] = refseq
    
    visualize_list.append(tmp)

    # 2.Fragment for 1st double crossover

    # 组装模板
    ref_screeningmarker_seq = uha + screeningmarker_seq + dha
    
    ## 2.1 uha
    tmp ={}
    tmp["name"] = "Fragment for 1st crossover" 
    tmp["detail"] = "PCR product. The upstream homologous arm (UHA), the downstream homologous arm (DHA), screening marker and primer-7/8 are visualized." 
    tmp["seq"] = ref_screeningmarker_seq
    tmp["params"] = [
        {
            'name':"uha",
            'start':str(ref_screeningmarker_seq.find(uha)),
            'end':str(ref_screeningmarker_seq.find(uha) + len(uha)),
            'direction':1,
            'color':'red',
        }
    ]
    
    ## 2.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(ref_screeningmarker_seq.find(dha)),
            "end":str(ref_screeningmarker_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )
    # 2.3 screeningmarker_seq
    tmp["params"].append(
        {
            "name":"screeningmarker",
            "start":str(ref_screeningmarker_seq.find(screeningmarker_seq)),
            "end":str(ref_screeningmarker_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
            'direction':1,
            'color':'pink',
        }
    )
    ## 2.3 primer7
    primer7 = dict_screeningmarker['PRIMER_LEFT_WHOLE_SEQUENCE']
    tmp["params"].append(
        {
            "name":"primer-7",
            "start":str(ref_screeningmarker_seq.find(primer7)),
            "end":str(ref_screeningmarker_seq.find(primer7) + len(primer7)),
            'direction':1,
            'color':'orange',
        }
    )

    
    ## 2.3 primer8
    primer8 = dict_screeningmarker['PRIMER_RIGHT_WHOLE_SEQUENCE']
    primer8_rev = revComp(primer8)
    tmp["params"].append(
        {
            "name":"primer-8",
            "start":str(ref_screeningmarker_seq.find(primer8_rev)),
            "end":str(ref_screeningmarker_seq.find(primer8_rev) + len(primer8_rev)),
            'direction':-1,
            'color':'orange',
        }
    )
    
    visualize_list.append(tmp)

    # 3.Fragment for 2nd double crossover
    # for alter sequence, manipulation sequence
    ref_manipulation_seq = uha + alt + dha
    
    ## 3.1 uha
    tmp ={}
    tmp["name"] = "Fragment for 2nd crossover" 
    tmp["detail"] = "PCR product or synthesized ssDNA. The upstream homologous arm (UHA), the downstream homologous arm (DHA), inserted sequence and primer-5/6 (for long fragment insertion or substitution or just one single ssDNA for short inserted sequence) are visualized." 
    tmp["seq"] = ref_manipulation_seq
    tmp["params"] = [
        {
            'name':"uha",
            'start':str(ref_manipulation_seq.find(uha)),
            'end':str(ref_manipulation_seq.find(uha) + len(uha)),
            'direction':1,
            'color':'red',
        }
    ]
    
    ## 3.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(ref_manipulation_seq.find(dha)),
            "end":str(ref_manipulation_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )

    # 3.3 seq_altered
    if mutation_type != "deletion":
        tmp["params"].append(
            {
                "name":"seq_altered",
                "start":str(ref_manipulation_seq.find(uha) + len(uha)),
                "end":str(ref_manipulation_seq.find(dha)),
                'direction':1,
                'color':'pink',
            }
        )
    # 增加ssDNA
    if dict_input_seq["level"] == "<=length_homologous_sequence":
        tmp["params"].append(
        {
            "name":"ssDNA",
            "start":str(ref_manipulation_seq.find(uha)),
            "end":str(ref_manipulation_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'#6c757d', 
        }
    )

    # 插入序列的引物
    if dict_insert:
        # 3.5 primer5
        primer5 = dict_insert['PRIMER_LEFT_WHOLE_SEQUENCE']
        tmp["params"].append(
            {
                "name":"primer-5",
                # "start":str(ref_manipulation_seq.find(primer5)),
                "start":str(ref_manipulation_seq.find(uha) + len(uha)),
                # "end":str(ref_manipulation_seq.find(primer5) + len(primer5)),
                "end":str(ref_manipulation_seq.find(uha) + len(uha) + len(primer5)),
                'direction':1,
                'color':'#0dcaf0',
            }
        )

        # 3.6 primer6
        primer6 = dict_insert['PRIMER_RIGHT_WHOLE_SEQUENCE']
        primer6_rev = revComp(primer6)
        tmp["params"].append(
            {
                "name":"primer-6",
                # "start":str(ref_manipulation_seq.find(primer6_rev)),
                "start":str(ref_manipulation_seq.find(dha) - len(primer6)),
                # "end":str(ref_manipulation_seq.find(primer6_rev) + len(primer6_rev)),
                "end":str(ref_manipulation_seq.find(dha)),
                'direction':-1,
                'color':'#0dcaf0',
            }
        )
    else:
        primer5 =""
        primer6 =""

    visualize_list.append(tmp)
    
    # 4. verify1 visualize
    # verify1 验证
    v1_seq = dict_input_seq["seq_uha_verify"]+ screeningmarker_seq + dict_input_seq["seq_dha_verify"]
    ## 4.1 uha
    tmp = {}
    tmp["name"]="Verify1"
    tmp["detail"]="Verification of the 1st-round of double crossover and isolation. The upstream homologous arm (UHA), the downstream homologous arm (DHA) screening marker and test-primer-1/2 are visualized."
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
    
    ## 4.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(v1_seq.find(dha)),
            "end":str(v1_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )
    ## 4.3 screeningmarker
    tmp["params"].append(
        {
            "name":"screeningmarker",
            "start":str(v1_seq.find(screeningmarker_seq)),
            "end":str(v1_seq.find(screeningmarker_seq) + len(screeningmarker_seq)),
            'direction':1,
            'color':'pink',
        }
    )
    ## 4.4 v1 primer1
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
    
    
    ## 4.5 v1 primer2
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
    
    # 5. verify2 visualize
    # verify2 验证
    v2_seq = dict_input_seq["seq_uha_verify"]+ alt + dict_input_seq["seq_dha_verify"]
    
    ## 5.1 uha
    tmp ={}
    tmp["name"] = "Verify2"
    tmp["detail"] = "Verification of the 2nd-round of double crossover and isolation. The upstream homologous arm (UHA), the downstream homologous arm (DHA) and test-primer-3/4 are visualized."
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
    
    ## 5.2 dha
    tmp["params"].append(
        {
            "name":"dha",
            "start":str(v2_seq.find(dha)),
            "end":str(v2_seq.find(dha) + len(dha)),
            'direction':1,
            'color':'red',
        }
    )
    ## 5.3 Inserted fragment
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
    ## 5.4 v2 primer1
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
    
    ## 5.5 v2 primer2
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
        primer5,
        primer6,
        primer7,
        primer8,
        dict_verify1['PRIMER_LEFT_0_SEQUENCE'],
        dict_verify1['PRIMER_RIGHT_0_SEQUENCE'],
        dict_verify2['PRIMER_LEFT_0_SEQUENCE'],
        dict_verify2['PRIMER_RIGHT_0_SEQUENCE'],
    ]
    fsave.write(
        '\t'.join(linelist) + "\n"
    )

#Design process
def design_process(input_file_path,screeningmarker_file_path,workdir,ref_genome,config,assigned_tp1_seq=None,screening_maker_removal=None):
    dict_input_seq=input_to_primer_template(input_file_path,ref_genome,config,workdir)
    if isinstance(dict_input_seq,str):
        #print('*********** error ***************')
        #print(dict_input_seq)
        return dict_input_seq
    dict_primers_whole={'successful':{'verify1':[],'verify2':[],'inserted':[],'screening_marker':[]},'failed':{'verify1':[],'verify2':[],'inserted':[],'screening_marker':[]}}
    screeningmarker_seq=linearSeq(screeningmarker_file_path)
    print(screeningmarker_seq)
    failed_sheet_title=['PRIMER_LEFT_EXPLAIN','PRIMER_RIGHT_EXPLAIN','PRIMER_PAIR_EXPLAIN','PRIMER_LEFT_NUM_RETURNED','PRIMER_RIGHT_NUM_RETURNED','PRIMER_INTERNAL_NUM_RETURNED','PRIMER_PAIR_NUM_RETURNED']
    # 读写文件
    blast_input_file_path=workdir+'/blast_input.txt'
    primer_order_file_path=workdir+'/Primer_order.txt'
    visualize_file_path=workdir+'/Visualize_order.tsv'
    blast_input_file =  open(blast_input_file_path,'w')
    primer_order_fsave =  open(primer_order_file_path,'w')
    visualize_fsave =  open(visualize_file_path,'w')
    # write visualize file header
    headers = ['seqid','ref','alt',"manipulation","primer-5","primer-6","primer-7","primer-8","test-primer-1","test-primer-2","test-primer-3","test-primer-4"]
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
            dict_screening_marker_primers_attribute = screeningmarkerPrimerAttribute(key1,dictscreeningmarkerprimers,dict_input_seq[key1],screeningmarker_seq)
            dict_primers_whole['successful']['screening_marker'].append(dict_screening_marker_primers_attribute)
            target_seq=dict_input_seq[key1]['seq_uha']+dict_input_seq[key1]['seq_altered']+dict_input_seq[key1]['seq_dha']
            homo_arm_seq_1=dict_input_seq[key1]['seq_uha']
            homo_arm_seq_2=dict_input_seq[key1]['seq_uha']
            blast_input_file.write('>'+key1+'_UHA'+'\n'+homo_arm_seq_1+'\n')
            blast_input_file.write('>'+key1+'_DHA'+'\n'+homo_arm_seq_2+'\n')
            if dict_input_seq[key1]["level"] == ">length_homologous_sequence":
                inserted_temp = dict_input_seq[key1]['seq_altered']
                dictinsertedprimers=primerDesign(key1,inserted_temp,config,"insert")
                if len(dictinsertedprimers)<10:
                    dict_inserted_failed_output={k:dictinsertedprimers[k] for k in dictinsertedprimers.keys() if k in failed_sheet_title}
                    dict_inserted_failed_output['ID']=key1
                    dict_primers_whole['failed']['inserted'].append(dict_inserted_failed_output)
                    print("error")
                else:
                    dict_inserted_primers_attribute = insertPrimerAttribute(key1,dictinsertedprimers,dict_input_seq[key1],inserted_temp)
                    dict_primers_whole['successful']['inserted'].append(dict_inserted_primers_attribute)
                    seq_verify1_temp = screeningmarker_seq[-400:] + dict_input_seq[key1]["seq_dha_verify"]
                    dictverify1primers=primerDesign(key1,seq_verify1_temp,config,"verify1")
                    if len(dictverify1primers)<10:
                        dict_verify1_failed_output={k:dictverify1primers[k] for k in dictverify1primers.keys() if k in failed_sheet_title}
                        dict_verify1_failed_output['ID']=key1
                        dict_primers_whole['failed']['verify1'].append(dict_verify1_failed_output)
                        print("error")
                    else:
                        if assigned_tp1_seq:
                            seq_verify1_temp_add_screeningmarker=screeningmarker_seq + dict_input_seq[key1]["seq_dha_verify"]
                            dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp,seq_verify1_temp_add_screeningmarker,assigned_tp1_seq)
                        else:
                            dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp)
                        dict_primers_whole['successful']['verify1'].append(dict_verify1_primers_attribute)
                        seq_verify2_temp = dict_input_seq[key1]["seq_uha_verify"] + dict_input_seq[key1]['seq_altered'] + dict_input_seq[key1]["seq_dha_verify"]
                        dictverify2primers=primerDesign(key1,seq_verify2_temp,config,"verify2")
                        if len(dictverify2primers)<10:
                            dict_verify2_failed_output={k:dictverify2primers[k] for k in dictverify2primers.keys() if k in failed_sheet_title}
                            dict_verify2_failed_output['ID']=key1
                            dict_primers_whole['failed']['verify2'].append(dict_verify2_failed_output)
                            print("error")
                        else:
                            dict_verify2_primers_attribute = verify2PrimerAttribute(key1,dictverify2primers,seq_verify2_temp)
                            dict_primers_whole['successful']['verify2'].append(dict_verify2_primers_attribute)
                            primers_submitted_output(dict_screening_marker_primers_attribute,dict_verify1_primers_attribute,primer_order_fsave,dict_verify2_primers_attribute,dict_inserted_primers_attribute)
            else:
                seq_verify1_temp = screeningmarker_seq[-400:] + dict_input_seq[key1]["seq_dha_verify"]
                dictverify1primers=primerDesign(key1,seq_verify1_temp,config,"verify1")
                if len(dictverify1primers)<10:
                    dict_verify1_failed_output={k:dictverify1primers[k] for k in dictverify1primers.keys() if k in failed_sheet_title}
                    dict_verify1_failed_output['ID']=key1
                    dict_primers_whole['failed']['verify1'].append(dict_verify1_failed_output)
                    print("error")
                else:
                    if assigned_tp1_seq:
                        seq_verify1_temp_add_screeningmarker=screeningmarker_seq + dict_input_seq[key1]["seq_dha_verify"]
                        dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp,seq_verify1_temp_add_screeningmarker,assigned_tp1_seq)
                    else:
                        dict_verify1_primers_attribute = verify1PrimerAttribute(key1,dictverify1primers,seq_verify1_temp)
                    dict_primers_whole['successful']['verify1'].append(dict_verify1_primers_attribute)
                    if dict_input_seq[key1]['type']=="deletion" and screening_maker_removal=="No":
                        primers_submitted_output(dict_screening_marker_primers_attribute,dict_verify1_primers_attribute,primer_order_fsave)
                    else:
                        seq_verify2_temp = dict_input_seq[key1]["seq_uha_verify"] + dict_input_seq[key1]['seq_altered'] + dict_input_seq[key1]["seq_dha_verify"]
                        dictverify2primers=primerDesign(key1,seq_verify2_temp,config,"verify2")
                        if len(dictverify2primers)<10:
                            dict_verify2_failed_output={k:dictverify2primers[k] for k in dictverify2primers.keys() if k in failed_sheet_title}
                            dict_verify2_failed_output['ID']=key1
                            dict_primers_whole['failed']['verify2'].append(dict_verify2_failed_output)
                            print("error")
                        else:
                            dict_verify2_primers_attribute = verify2PrimerAttribute(key1,dictverify2primers,seq_verify2_temp)
                            dict_primers_whole['successful']['verify2'].append(dict_verify2_primers_attribute)
                            primers_submitted_output(dict_screening_marker_primers_attribute,dict_verify1_primers_attribute,primer_order_fsave,dict_verify2_primers_attribute)
            # visualize
            '''
            generate_visualize_file(
                key1,
                dict_input_seq[key1],
                screeningmarker_seq,
                config,
                dict_screening_marker_primers_attribute,
                dict_verify1_primers_attribute,
                dict_verify2_primers_attribute,
                visualize_fsave,
                workdir,
                dict_inserted_primers_attribute,
                )
            '''
    #print(dict_primers_whole['successful']['dha'])
    # 保存文件
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
    parser.add_argument('--screeningmarker', '-s', help='screening marker seq file', required=True)
    parser.add_argument('--ref', '-r', help='reference genome', required=True)
    parser.add_argument('--dir', '-d', help='outputdir', required=True)
    parser.add_argument('--conf', '-c', help='config', required=True)
    parser.add_argument('--ustp', '-u', help='user specific test primer', required=False)
    parser.add_argument('--smr', '-m', help='screening maker removal', required=False)
    arguments = parser.parse_args()
    input_file_path = arguments.input
    screeningmarker_file_path = arguments.screeningmarker
    ref_genome = arguments.ref
    workdir = arguments.dir
    point_dict = conf_read(arguments.conf)
    user_specific_test_primer = arguments.ustp
    screening_maker_removal = arguments.smr
    time1=time.time()
    design_process(input_file_path,screeningmarker_file_path,workdir,ref_genome,point_dict,user_specific_test_primer,screening_maker_removal)
    time2=time.time()
    print(time2-time1)