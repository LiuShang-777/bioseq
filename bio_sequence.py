# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 11:20:44 2019

@author: liushang
"""
print('>==Starting the biosequence==<')
import argparse
parser=argparse.ArgumentParser(description='the option for fasta sequence')
parser.add_argument('-f','--fasta_file',help='the name of the input fasta file',type=str)
parser.add_argument('-n','--name_file',help='the name file in txt format',type=str)
parser.add_argument('-rf','--result_file',type=str,help='the directory of result file')
parser.add_argument('-c','--class_type',type=str,help='the class choose to opt for fasta',choices=['reverse','dna_to_rna','rna_to_dna','reverse_comp'])
parser.add_argument('-in','--inquery',type=str,help='the inquery fasta file')
parser.add_argument('-op','--option',type=str,help='the options you choose',choices=['blast','base','get_seq','find_motif','translate'])
parser.add_argument('-mtf','--motif',type=str,help='the motif you need to find in fasta file')
args=parser.parse_args()
inputfile_fasta=args.fasta_file
result_file=args.result_file
with open(inputfile_fasta,'r') as file:
    name=[]
    seq=[]
    sub_seq=[]
    for line in file:
        line=line.strip()
        if line[0]=='>':
            name.append(line[1:])
            seq.append(''.join(sub_seq))
            sub_seq=[]
        else:
            sub_seq.append(line)
    seq.append(''.join(sub_seq))
    seq=seq[1:]
    if len(name)!=len(seq):
        exit('the fasta file\'s format is incorrected')
    else:
        print('the input fasta file has been correctly imported')
if args.option=='base':    
    def basic_info(result_file):
        dic={}
        for i,j in zip(name,seq):
            dic[i]=[]
            length=len(j)
            cg_contnet=int((j.count('C')+j.count('G'))*100/length)
            atg_count1,atg_count2,atg_count3=0,0,0
            atg_loc1,atg_loc2,atg_loc3=[],[],[]
            for ele in range(0,len(j),3):
                sub_str=j[ele:ele+3]
                if (sub_str=='ATG')|(sub_str=='AUG'):
                    atg_count1+=1
                    atg_loc1.append(ele)
            for ele in range(1,len(j),3):
                sub_str=j[ele:ele+3]
                if (sub_str=='ATG')|(sub_str=='AUG'):
                    atg_count2+=1
                    atg_loc2.append(ele)
            for ele in range(2,len(j),3):
                sub_str=j[ele:ele+3]
                if (sub_str=='ATG')|(sub_str=='AUG'):
                    atg_count3+=1
                    atg_loc3.append(ele)
            dic[i]=[length,cg_contnet,atg_count1,atg_loc1,atg_count2,atg_loc2,atg_count3,atg_loc3]
        with open(result_file,'w') as file:
            for i in dic.keys():     
                file.write('gene: %s\n'%i)
                file.write('length: %d\n'%dic[i][0])
                file.write('GC content: %d\n'%dic[i][1])
                file.write('number of ATG in s1 : %d\n'%dic[i][2])
                file.write('location of ATG in s1 :\n')
                file.writelines(str(dic[i][3])+'\n')
                file.write('number of ATG in s2 : %d\n'%dic[i][4])
                file.write('location of ATG in s2 :\n')
                file.writelines(str(dic[i][5])+'\n')
                file.write('number of ATG in s3 : %d\n'%dic[i][6])
                file.write('location of ATG in s3 :\n')
                file.writelines(str(dic[i][7])+'\n')
    basic_info(result_file+'basicinfo.txt')  
elif args.option=='find_motif':
    motif=args.motif
    def find_motif(motif,result_file) :
        dic={}
        for i,j in zip(name,seq):
            dic[i]=[]
            motif_count1,motif_count2,motif_count3=0,0,0
            motif_loc1,motif_loc2,motif_loc3=[],[],[]
            for ele in range(0,len(j),3):
                sub_str=j[ele:ele+3]
                if sub_str==motif:
                   motif_count1+=1
                   motif_loc1.append(ele)
            for ele in range(1,len(j),3):
                sub_str=j[ele:ele+3]
                if sub_str==motif:
                   motif_count2+=1
                   motif_loc2.append(ele)
            for ele in range(2,len(j),3):
                sub_str=j[ele:ele+3]
                if sub_str==motif:
                   motif_count3+=1
                   motif_loc3.append(ele)
            dic[i]=[motif_count1,motif_loc1,motif_count2,motif_loc2,motif_count3,motif_loc3]
        with open(result_file,'w') as file:
            for i in dic.keys():
                file.write('gene: %s\n'%i)
                file.write('motif_count in s1: %d\n'%dic[i][0])
                file.write('motif_loc in s1:\n')
                file.writelines(str(dic[i][1])+'\n')
                file.write('motif_count in s2: %d\n'%dic[i][2])
                file.write('motif_loc in s2:\n')
                file.writelines(str(dic[i][3])+'\n')
                file.write('motif_count in s3: %d\n'%dic[i][4])
                file.write('motif_loc in s3:\n')
                file.writelines(str(dic[i][5])+'\n')
    find_motif(motif,result_file+'motif_find.txt')
elif args.option=='get_seq':
    inputfile_txt= args.name_file
    def get_sequence(inputfile_txt,result_file):
        global name,seq
        with open(inputfile_txt,'r') as file:
            gene_name=[]
            get_seq=[]
            get_name=[]
            for line in file:
                line=line.strip()
                gene_name.append(line)
            for i in gene_name:
                if i in name:
                    get_name.append(i)
                    get_seq.append(seq[name.index(i)])
                else:
                    print('%s is not in fasta file'%i)
        with open(result_file,'w') as file:
            for i,j in zip(get_name,get_seq):
                file.write('>'+i+'\n')
                file.write(j+'\n')
    get_sequence(inputfile_txt,result_file+'result_get_seq.fa')
elif args.option=='blast':
    inquery=args.inquery
    inquery_name=[]
    inquery_seq=[]
    inquery_sub_seq=[]
    with open(inquery,'r') as file:
        for line in file:
            line=line.strip()
            if line[0]=='>':
                inquery_name.append(line)
                inquery_seq.append(''.join(inquery_sub_seq))
                inquery_sub_seq=[]
            else:
                inquery_sub_seq.append(line)
        inquery_seq.append(''.join(inquery_sub_seq))
        inquery_seq=inquery_seq[1:]
        if len(inquery_name)==len(inquery_seq):
            print('inquery file has been imported')
        else:
            print('the format of inquery file is wrong')
    def letter_score(a,b):
        two_str=a+b
        if a==b:
            return 2
        elif (('A' in two_str)&('G' in two_str))|(('C' in two_str)&('T' in two_str)):
            return -5
        else:
            return -7
    def blast(sequence_in,sequence_db):
        gap=-5
        score_seq_in=[]
        point_seq_in=[]
        for i in range(len(sequence_in)+1):
            if i==0:
                sub_point_seq_in,sub_score_seq_in=[0],[0]
                for i in range(1,len(sequence_db)+1):
                    sub_score_seq_in.append(gap*1)
                    sub_point_seq_in.append('db_gap')
            else:
                sub_score_seq_in,sub_point_seq_in=[],[]
                sub_score_seq_in.append(gap*i)
                sub_point_seq_in.append('in_gap')
            score_seq_in.append(sub_score_seq_in)
            point_seq_in.append(sub_point_seq_in)
        for i in range(1,len(sequence_in)+1):
            nucl1=sequence_in[i-1]
            for j in range(1,len(sequence_db)+1):
                nucl2=sequence_db[j-1]
                score=letter_score(nucl1,nucl2)
                dia_score=score_seq_in[i-1][j-1]+score
                db_gap_score=score_seq_in[i][j-1]+gap
                in_gap_score=score_seq_in[i-1][j]+gap
                select_score=max(dia_score,db_gap_score,in_gap_score)
                score_seq_in[i].append(select_score)
                if score_seq_in[i][j]==dia_score:
                    point_seq_in[i].append('normal')
                elif score_seq_in[i][j]==db_gap_score:
                    point_seq_in[i].append('db_gap')
                else:
                    point_seq_in[i].append('in_gap')
        seq_in=''
        seq_db=''
        i,j=len(sequence_in),len(sequence_db)
        while True:
            if point_seq_in[i][j]==0:
                break
            elif point_seq_in[i][j]=='normal':
                seq_in+=sequence_in[i-1]
                seq_db+=sequence_db[j-1]
                i-=1
                j-=1
            elif point_seq_in[i][j]=='in_gap':
                seq_in+=sequence_in[i-1]
                seq_db+='_'
                i-=1
            else:
                seq_in+='_'
                seq_db+=sequence_db[j-1]
                j-=1
        seq_in=seq_in[::-1]
        seq_db=seq_db[::-1]
        print(seq_in+'\n')
        print(seq_db+'\n')
    for i,j,a,b in zip(seq,inquery_seq,name,inquery_name):
        print('%s vs %s\n'%(a,b))
        blast(i,j)
elif args.option=='translate':
    class_type=args.class_type                
    def translate(class_type,result_file):
        result_seq=[]
        if class_type=='reverse':
            for i in seq:
                i=i[::-1]
                result_seq.append(i)
        elif class_type=='dna_to_rna':
            for i in seq:
                i=i.replace('T','U')
                result_seq.append(i)
        elif class_type=='rna_to_dna':
            for i in seq:
                i=i.replace('U','T')
                result_seq.append(i)
        elif class_type=='reverse_comp':
            trans=str.maketrans('ATGC','TACG')
            for i in seq:
                i=i.translate(trans)
                i=i[::-1]
                result_seq.append(i)
        with open(result_file,'w') as file:        
            for i,j in zip(name,result_seq):
                file.write('>'+i+'\n')
                file.write(j+'\n')
    translate(class_type,result_file+'translate.txt')
print('the biosequence analysis has finished')               