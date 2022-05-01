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

def compl_seq(strseq):
    strseq=strseq.upper()
    c_strseq=strseq.translate(string.maketrans("ATCG", "TAGC"))
    return c_strseq

def find_head_index(myref,substr):
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
    for i in range(sLen-pLen):
        tmp_str=myref[i:i+pLen]
        MM=Levenshtein.hamming(substr, tmp_str)
        hm.append([i,MM])
    if len(hm)>0:
        hm=sorted(hm, key=lambda hm:hm[1])
        if hm[0][1] <= maxMM:
            return_pos=hm[0][0] 
    return return_pos

def find_tail_index(myref,substr1):
    myref=myref[::-1]
    substr=substr1[::-1]
    myindex=find_head_index(myref, substr)
    if myindex>=0:
        return([myindex])
    elif len(substr)>7 and substr[-7:] in myref:
        temp_index=myref.index(substr[-7:])
        if temp_index>=0:
            match_len=7
            while temp_index>0:
                newlen=match_len+temp_index
                if newlen<=len(substr) and newlen<=len(myref):
                    temp_substr=substr[-newlen:]
                    temp_index1=find_head_index(myref, temp_substr)
                    if temp_index1>=0:
                        match_len=len(temp_substr)
                        temp_index=temp_index1
                    else:
                        break
                else:
                    break
            if temp_index<3:
                if match_len*2<len(substr):
                    sub1=substr1[:match_len]
                    sub2=substr1[match_len:]
                    n=find_head_index(sub2, sub1)
                    if n>=0:
                        match_len1=match_len*2+n
                        temp_index=find_head_index(myref, substr[-match_len1:])
                        if temp_index>=0:
                            match_len=match_len1
                            while temp_index>0:
                                newlen=match_len+temp_index
                                if newlen<=len(substr) and newlen<=len(myref):
                                    temp_substr=substr[-newlen:]
                                    temp_index1=find_head_index(myref, temp_substr)
                                    if temp_index1>=0:
                                        match_len=len(temp_substr)
                                        temp_index=temp_index1
                                    else:
                                        break
                                else:
                                    break
                
                return([len(substr),match_len])
            else:
                return([-1])
    else:
        return([-1])

def find_STR_F(temp_str):
    STR_len=0
    temp_str_len=len(temp_str)
    outlist=[]
    out1=[]
    line=temp_str
    for pname in pnames:
        FP,RP,FP_len,RP_len,motif,motif_len,ref,ref_len=p_info_dict[pname]
        FP1,RP1=FP,RP
        m=find_head_index(temp_str, FP)
        if RP:
            RP_len=int(RP_len)
        else:
            RP_len=0
        if FP:
            FP_len=int(FP_len)
        else:
            FP_len=0
        if motif:
            motif=motif.split("/")[0]
        if m>=0:
            #print("cunzai FP")
            FP1="FP"
            n=find_tail_index(temp_str, RP)
            if len(n)>1:
                RP1="partial"
                if RP_len==n[0]:
                    #print("match partial RP")
                    STR_len=temp_str_len-m+(n[0]-n[1])
                    sub_str=temp_str[m+FP_len:-n[1]]
                    if motif*2 in sub_str:
                        outlist.append("%s\t%s\t%s\t%s\t%s"%(pname,str(STR_len),str(len(temp_str)),FP+sub_str+RP,"\t".join(p_info_dict[pname])))
                else: 
                    print("RP_len not eq: %s\n%s\t%s\n"%(line,pname,"\t".join(p_info_dict[pname])))
                    continue
            elif len(n)==1:
                if n[0]>=0:
                    RP1="RP"
                    #print("cunzai RP")
                    STR_len=temp_str_len-m-n[0]
                    sub_str=temp_str[m+FP_len:-n[0]-RP_len]
                    if motif*2 in sub_str:
                        outlist.append("%s\t%s\t%s\t%s\t%s"%(pname,str(STR_len),str(len(temp_str)),FP+sub_str+RP,"\t".join(p_info_dict[pname])))
                    
            out1.append("%s\t%s\t%s\t%s"%(pname,FP1,RP1,STR_len))
    #print(outlist)
    #print("match primer: %s"%len(outlist))
    #print(out1)
    if len(outlist)==1:
        return(outlist,out1)
    elif len(outlist)>1:
        return(outlist,out1)
    else:
        return(-1,out1)

def find_STR_R(temp_str):
    STR_len=0
    line_len=len(temp_str)
    outlist=[]
    out1=[]
    line=temp_str
    for pname in pnames:
        FP,RP,FP_len,RP_len,motif,motif_len,ref,ref_len=p_info_dict[pname]
        FP1,RP1=FP,RP
        if RP:
            RP_len=int(RP_len)
            RP=reverse_compl_seq(RP)
        else:
            RP_len=0
        if FP:
            FP_len=int(FP_len)
            FP=reverse_compl_seq(FP)
        else:
            FP_len=0
        if motif:
            motif_len=int(motif_len)
            motif=reverse_compl_seq(motif)
        else:
            motif_len=0
        m=find_head_index(line, RP)
        if m>=0:
            #print("cunzai RP_rc")
            RP1="RP"
            n=find_tail_index(line, FP)
            if len(n)>1:
                FP1="partial"
                if FP_len==n[0]:
                    #print("match partial FP")
                    STR_len=line_len-m+(n[0]-n[1])
                    sub_str=line[m+RP_len:-n[1]]
                    if motif*2 in sub_str:
                        outlist.append("%s\t%s\t%s\t%s\t%s"%(pname,str(STR_len),str(len(line)),reverse_compl_seq(RP+sub_str+FP),"\t".join(p_info_dict[pname])))
                    
            elif len(n)==1:
                if n[0]>=0:
                    FP1="FP"
                    #print("cunzai FP_rc")
                    STR_len=line_len-m-n[0]
                    sub_str=line[m+RP_len:-n[0]-FP_len]
                    if motif*2 in sub_str:
                        outlist.append("%s\t%s\t%s\t%s\t%s"%(pname,str(STR_len),str(len(line)),reverse_compl_seq(RP+sub_str+FP),"\t".join(p_info_dict[pname])))
                    
            out1.append("%s\t%s\t%s\t%s"%(pname,FP1,RP1,STR_len))
    #print("match primer: %s"%len(outlist))
    #print(outlist)
    #print(out1)
    if len(outlist)==1:
        return(outlist,out1)
    elif len(outlist)>1:
        return(outlist,out1)
    else:
        return(-1,out1)
    
    
def find_STR(line):
    if line!="":
        num, temp_str=line.split()
        temp_str=temp_str.upper()
        num=int(num)
        strlen=len(temp_str)
        #print(line)
        F_str,fout1=find_STR_F(temp_str)
        if F_str>0:
            #print(F_str+[num])
            return(F_str+[num])
        else:
            R_str,rout1=find_STR_R(temp_str)
            if R_str>0:
                return(R_str+[num])
            else:
                return(-1)

      
pfile=open(options.info_file,'r')
plines=pfile.read().split("\n")
p_info_dict={}
for pline in plines[1:-1]:
    infos=pline.split("\t")
    p_info_dict[infos[0]]=infos[1:]

pnames=p_info_dict.keys()
print(type(pnames),len(pnames))
pnames=sorted(pnames, reverse=True)
mult=0
outfile=open(options.outfile,"w")
fqfile=open(options.reads_file,"r")
lines=fqfile.read().split("\n")

AMELX="ATAGTGTGTTGATTCTTTATCCCAGATGTTTCTCAAGTGGTCCTGATTTTACAGTTCCTACCACCAGCTTC"
AMELY="ATAGTGGGTGGATTCTTCATCCCAAATAAAGTGGTTTCTCAAGTGGTCCCAATTTTACAGTTCCTACCATCAGCTTC"
Yindel_1="GATTTAAACTCTCTGAATCAGGCACATGCCTTCTCACTTCTCAAGAATGAACAG"
Yindel_2="GATTTAAACTCTCTGAATCAGGCACATGCCTTCTCACTTCTCTTCTCAAGAATGAACAG"
def spl_site(tmpstr,site):
    tmpstr=tmpstr.upper()
    tmpstr_rc=reverse_compl_seq(tmpstr)
    strlen=len(tmpstr)
    a=find_head_index(tmpstr, site)
    if a>=0:
        return(a)
    else:
        b=find_head_index(tmpstr_rc, site)
        if b>=0:
            return(b)
        else:
            return(-1)

if __name__=="__main__":
    print("reads num: %s"%len(lines))
    amx=0
    amy=0
    y1=0
    y2=0
    pool=Pool(processes=thread_number)
    find_STR_results=[]
    for seq in lines[:-1]:
        num1, temp_str=seq.split()
        num1=int(num1)
        a1=spl_site(temp_str, AMELX)
        if a1>0:
            amx=amx+num1
        a2=spl_site(temp_str, AMELY)
        if a2>0:
            amy=amy+num1 
        a3=spl_site(temp_str, Yindel_1)
        if a3>0:
            y1=y1+num1 
        a4=spl_site(temp_str, Yindel_2)
        if a4>0:
            y2=y2+num1
        find_STR_results.append(pool.apply_async(find_STR, args=(seq,)))
    pool.close()
    pool.join()
    site_len={}
    str_reads=[]
    for a in find_STR_results:
        reg1=a.get()
        if reg1>0:
            if len(reg1)>2:
                print("match mult site: ",reg)
            elif len(reg1)==2:
                reg=reg1[0].split("\t")
                name=reg[0]
                s_len=int(reg[1])
                site_len.setdefault(name, []).append([s_len,int(reg1[-1]),reg[3]])
                #outfile.write("%s\t%s\t%s\t%s\t%s\n"%(reg[0],reg[1],reg[2],str(reg1[-1]),"\t".join(reg[3:])))
    print("site num: %s"%len(site_len))
    outfile.write("AMEL\t\tX\t%s\t%s\nAMEL\t\tY\t%s\t%s\nYindel\t\t1\t%s\t%s\nYindel\t\t2\t%s\t%s\n"%(amx,AMELX,amy,AMELY,y1,Yindel_1,y2,Yindel_2))
    sites=site_len.keys()
    site_genotype={}
    for site in sites:
        temp_list=site_len[site]
        len_depth={}
        len_str={}
        str_depth={}
        for temp in temp_list:
            if int(temp[0])==len(temp[2]):
                str_depth.setdefault(temp[2],[]).append(temp[1])        
                len_depth.setdefault(temp[0],[]).append(temp[1])
            else:
                print(site,temp[0],len(temp[2]),temp[2])
        lenkeys=len_depth.keys()
        for s in str_depth.keys():
            if len(s) in lenkeys:
                len_str.setdefault(len(s),[]).append([sum(str_depth[s]),s])
            else:
                len_str[len(s)]=[sum(str_depth[s]),s]
        d=sorted(len_depth.items(), key=lambda len_depth:len_depth[0])
        gt={}
        #print(site)
        depth=[]
        if p_info_dict[site][4]:
            motif_len=int(p_info_dict[site][5])
            alt=int(p_info_dict[site][7])
            alt_len=int(p_info_dict[site][6])
            for k in d:
                #print(site,k,len_str)
                temp_gt=""
                if k[0]==alt_len:
                    temp_gt=alt
                elif k[0]>alt_len:
                    temp_gt1=alt+(k[0]-alt_len)/motif_len
                    temp_gt2=(k[0]-alt_len)%motif_len
                    if temp_gt2:
                        temp_gt="%s.%s"%(temp_gt1,temp_gt2)
                    else:
                        temp_gt=temp_gt1
                elif k[0]<alt_len:
                    temp_gt1=alt-int(math.ceil((alt_len-k[0])*1.0/motif_len))
                    temp_gt2=(alt_len-k[0])%motif_len 
                    if temp_gt2:
                        temp_gt="%s.%s"%(temp_gt1,motif_len-temp_gt2)
                    else:
                        temp_gt=temp_gt1   
                gt[temp_gt]=sum(k[1])
                depth.append(sum(k[1]))
                if len(len_str[k[0]])==1:
                    outfile.write("%s\t%s\t%s\t%s\t%s:%s\n"%(site,k[0],temp_gt,sum(k[1]),len_str[k[0]][0][0],len_str[k[0]][0][1]))
                else:
                    tmp=len_str[k[0]]
                    tmp=sorted(tmp, key=lambda tmp:tmp[0], reverse=True)
                    outfile.write("%s\t%s\t%s\t%s\t"%(site,k[0],temp_gt,sum(k[1])))
                    for tmp1 in tmp:
                        outfile.write("%s:%s;"%(tmp1[0],tmp1[1]))
                    outfile.write("\n")
        else:
            for k in d:
                print(site,k)
                if len(len_str[k[0]])==1:
                    outfile.write("%s\t%s\t%s\t%s\t%s\n"%(site,k[0],temp_gt,sum(k[1]),len_str[k[0]]))
                else:
                    outfile.write("%s\t%s\t\t%s\t%s\n"%(site,k[0],sum(k[1]),";".join(len_str[k[0]])))
                depth.append(sum(k[1]))
        outfile.write("%s\ttotal depth\t\t%s\t\n"%(site,sum(depth)))

outfile.close()
fqfile.close()
pfile.close()
#outfile1.close()
endtime=datetime.datetime.now()
print(time.strftime("%y/%m/%d %H:%M:%S",time.localtime(time.time())))
print("time used :",str(endtime-starttime))
print("well done!")    
