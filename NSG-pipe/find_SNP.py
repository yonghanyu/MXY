from multiprocessing import Process as pro
from multiprocessing import Pool
from optparse import OptionParser
from collections import Counter
import string,time,datetime
import Levenshtein
import os,math

starttime=datetime.datetime.now()
print(time.strftime("%y/%m/%d %H:%M:%S",time.localtime(time.time())))

parser=OptionParser()
parser.set_defaults(thread_number=2)
parser.add_option("-i","--info_file", dest="info_file", help="Input your primer info file in FP,RP,FF,RF format", metavar="PRIMER")
parser.add_option("-f","--reads_file", dest="reads_file", help="Input uniq reads file in decompression format", metavar="reads_file")
parser.add_option("-o","--out_file", dest="outfile", help="Input outfile name", metavar="outfile")
parser.add_option("-p","--threads", dest="thread_number", help="Input the thread number(default:2)", metavar="THREADS")
(options,args)=parser.parse_args()

thread_number=int(options.thread_number)

def reverse_compl_seq(strseq):
    strseq=strseq.upper()
    rc_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))[::-1]
    return rc_strseq

def find_index(myref,substr):
    myindex=-1
    if substr in myref:
        myindex=myref.index(substr)
    else:
        myindex=FUZZYMATCH(myref, substr)
    return(myindex)

def FUZZYMATCH(myref,substr):
    MM=100
    startPos=-1
    pLen = len(substr)
    sLen = len(myref)
    hm=[]
    maxMM=pLen/15
    return_pos=-1
    for i in range(sLen-pLen+1):
        tmp_str=myref[i:i+pLen]
        MM=Levenshtein.hamming(substr, tmp_str)
        hm.append([i,MM])
    if len(hm)>0:
        hm=sorted(hm, key=lambda hm:hm[1])
        if hm[0][1] <= maxMM:
            return_pos=hm[0][0] 
    return return_pos

def find_SNP_F(temp_str):
    b=-1
    line_len=len(temp_str)
    outlist=[]
    for pname in pnames:
        sub_str,flen,rlen=p_info_dict[pname]
        a=find_index(temp_str, sub_str)
        if a>=0:
            b=a+int(flen)
            outlist.append([pname,a,temp_str[b],temp_str])
    if len(outlist)>1:
        print(len(outlist))
    return([b, outlist])

def find_SNP_R(temp_str):
    b=-1
    line_len=len(temp_str)
    outlist=[]
    line=reverse_compl_seq(temp_str)
    for pname in pnames:
        sub_str,flen,rlen=p_info_dict[pname]
        a=find_index(line, sub_str)
        if a>=0:
            b=a+int(flen)
            outlist.append([pname,a,line[b],line])
    if len(outlist)>0:
        print(len(outlist),b,outlist)
    return([b, outlist])    

def find_SNP(line):
    if line!="":
        num, temp_str=line.split()
        temp_str=temp_str.upper()
        num=int(num)
        strlen=len(temp_str)
        #print(line)
        F_str, site_list=find_SNP_F(temp_str)
        if F_str>0:
            #print(len(site_list))
            return([num]+site_list)
        else:
            R_str, site_list=find_SNP_R(temp_str)
            if R_str>0:
                #print(len(site_list))
                return([num]+site_list)
            else:
                return(-1)

      
pfile=open(options.info_file,'r')
#pfile=open(r"E:\FGI\data\20170808\STR\data\snp.near25.info.txt","r")
plines=pfile.read().split("\n")
p_info_dict={}
for pline in plines[0:-1]:
    infos=pline.split("\t")
    p_info_dict[infos[0]]=infos[1:]

pnames=p_info_dict.keys()

#fqfile=open(r"E:\FGI\data\20170808\STR\F2_R1_d10.reads","r")
fqfile=open(options.reads_file,"r")
lines=fqfile.read().split("\n")
outfile=open(options.outfile,"w")
outfile.write("%s Site\talt1\talt2\tdepth\tall alt\n"%os.path.basename(options.reads_file).split("_")[0])
if __name__=="__main__":
    print("reads num: %s"%len(lines))
    amx=0
    amy=0
    y=0
    pool=Pool(processes=thread_number)
    snp_depth={}
    find_SNP_results=[]
    for line in lines[:-1]:
        #out=find_SNP(line)
        find_SNP_results.append(pool.apply_async(find_SNP, args=(line,)))
    pool.close()
    pool.join()
    site_len={}
    str_reads=[]
    for a in find_SNP_results:
        reg1=a.get()
        if reg1>0:
            if len(reg1)==2:
                snp_depth.setdefault(reg1[1][0],[]).append([reg1[1][2],int(reg1[0])])
            else:
                print(reg1)
    print("site num: %s"%len(snp_depth))
                
    snps=snp_depth.keys()
    snp_gt={}
    for snp in snps:
        temp_list=snp_depth[snp]
        tmp_dict={}
        for temp in temp_list:
            tmp_dict.setdefault(temp[0],[]).append(temp[1])
        for k in tmp_dict.keys():
            tmp=tmp_dict[k]
            tmp_dict[k]=sum(tmp)
        d=sorted(tmp_dict.items(), key=lambda tmp_dict:tmp_dict[1], reverse=True)
        t=d[0][1]
        temp_gt1=d[0][0]
        temp_gt2=""
        outstr=["%s:%s"%(d[0][0],d[0][1]),]
        rate=[]
        alt=[]
        for d1 in d[1:]:
            t=t+d1[1]
            r=d1[1]*100.0/d[0][1]
            rate.append("%.2f"%r)
            outstr.append("%s:%s"%(d1[0],d1[1]))
            if r>15:
                alt.append(d1[0])
        if len(alt)>0:
            temp_gt2=d[1][0]
        elif len(rate)>0:
            if rate[0]>10 and rate[0]<15:
                temp_gt2="mb%s"%d[1][0]
        outfile.write("%s\t%s\t%s\t%s\t%s\n"%(snp,temp_gt1,temp_gt2,t,";".join(outstr)))
   
fqfile.close()
pfile.close()
outfile.close()
endtime=datetime.datetime.now()
print(time.strftime("%y/%m/%d %H:%M:%S",time.localtime(time.time())))
print("time used :",str(endtime-starttime))
print("well done!")    
