from optparse import OptionParser
from collections import Counter
import string,time,datetime
import os,math

starttime=datetime.datetime.now()
print(time.strftime("%y/%m/%d %H:%M:%S",time.localtime(time.time())))

parser=OptionParser()
parser.add_option("-i","--in_file", dest="in_file", help="Input your primer info file in FP,RP,FF,RF format", metavar="PRIMER")
parser.add_option("-o","--out_file", dest="outfile", help="Input outfile name", metavar="outfile")
(options,args)=parser.parse_args()

infile=open(options.in_file,"r")
outfile=open(options.outfile,"w")
gt_dict={}

lines=infile.read().split("\n")
for line in lines:
    items=line.split("\t")
    if len(items)==5 and "total depth" not in items[1]:
        if "AMEL" not in items[0]:
            if items[0]=="":
                if "TAGATAGATGATAGATAGA" in items[4].split(";")[0]:
                    gt_dict.setdefault(items[0],[]).append([float(items[2]),int(items[3]),items[4].split(";")[0]])
                elif int(items[3])>35:
                    gt_dict.setdefault(items[0],[]).append([float(items[2]),int(items[3]),items[4].split(";")[0]])
            else:
                gt_dict.setdefault(items[0],[]).append([float(items[2]),int(items[3]),items[4].split(";")[0]])
        else:
            gt_dict.setdefault(items[0],[]).append([items[2],int(items[3]),items[4].split(";")[0]])

gts=gt_dict.keys()
for gt in gts:
    site_gts=gt_dict[gt]
    site_gts=sorted(site_gts, key=lambda site_gts:site_gts[1], reverse=True)
    #depth=site_gts[0][1]
    if site_gts[0][1]>0:
        temp_gt1=site_gts[0][0]
        gtstr1=site_gts[0][2]
    else:
        temp_gt1=""
        gtstr1=""
    temp_gt2=""

    gtstr2=""
    if len(site_gts)==1:
        outfile.write("%s\t%s\t\t\t%s\t%s\n"%(gt,temp_gt1,site_gts[0][1],site_gts[0][2]))    
    elif len(site_gts)>1:
        if "AMEL" in gt or "Yindel" in gt:
            d1=""
            d2=""
            if int(site_gts[0][1])>0:
                temp_gt1=site_gts[0][0]
                d1=site_gts[0][1]
                gtstr1=site_gts[0][2]
            if int(site_gts[1][1])>0:
                d2=site_gts[1][1]
                temp_gt2=site_gts[1][0]
                gtstr2=site_gts[1][2]
            if temp_gt1 and temp_gt2:
                outfile.write("%s\t%s\t%s\t\t%s:%s\t%s;%s\n"%(gt,temp_gt1,temp_gt2,d1,d2,gtstr1,gtstr2))
            elif temp_gt1:
                outfile.write("%s\t%s\t\t\t%s\t%s\n"%(gt,temp_gt1,d1,gtstr1))
            elif temp_gt2:
                outfile.write("%s\t%s\t\t\t%s\t%s\n"%(gt,temp_gt2,d2,gtstr2))
        else: 
            temp=site_gts[0][1]
            rate=[]
            depthlist=["%s:%s"%(site_gts[0][0],site_gts[0][1]),]
            temp_gt1=site_gts[0][0]
            gtstr1=site_gts[0][2]
            temp_gt2=""
            alt=[]
            for site in site_gts[1:7]:
                rate.append("%.2f"%(site[1]*100.0/temp))
                depthlist.append("%s:%s"%(site[0],site[1]))
            
            if "DYS" not in gt:
                for site_gt in site_gts[1:]:
                    if site_gt[1]>0.16*temp:
                        alt.append(site_gt)
            elif  "ab" in gt:
                for site_gt in site_gts[1:]:
                    if site_gt[1]>0.16*temp:
                        alt.append(site_gt)
    
            if len(alt)==1:
                if float(rate[0])>30:
                    temp_gt2=alt[0][0]
                    gtstr2=alt[0][2]
                elif site_gts[1][1]>10:
                    if abs(float(temp_gt1)-float(alt[0][0]))>1:
                        temp_gt2=alt[0][0]
                        gtstr2=alt[0][2]
            elif len(alt)>1:
                if float(rate[0])>30:
                    if float(rate[1])<20 or float(alt[0][0])>float(temp_gt1):
                        temp_gt2=alt[0][0]
                        gtstr2=alt[0][2]
                elif site_gts[1][1]>10 and abs(float(temp_gt1)-float(alt[0][0]))>1:
                    temp_gt2=alt[0][0]
                    gtstr2=alt[0][2]
                elif site_gts[2][1]>10 and abs(float(temp_gt1)-float(alt[1][0]))>1:
                    temp_gt2=alt[1][0]
                    gtstr2=alt[1][2]
            elif len(alt)==0 :
                tmp1=[]
                if "DYS" not in gt or "ab" in gt:
                    tmp=sorted(site_gts, key=lambda site_gts:site_gts[0], reverse=True)
                    if tmp[0][1]>tmp[1][1]:
                        tmp1.append(tmp[0])
                    elif tmp[0][0]-tmp[1][0]>1:
                        tmp1.append(tmp[0])
                    for i in range(1,len(tmp)):
                        if abs(tmp[i][0]-tmp[i-1][0])==1:
                            if tmp[i][1]>tmp[i-1][1]:
                                if tmp[i] not in tmp1:
                                    tmp1.append(tmp[i])
                        elif abs(tmp[i][0]-tmp[i-1][0])>1:
                            if tmp[i] not in tmp1:
                                tmp1.append(tmp[i])
                        elif abs(tmp[i][0]-tmp[i-1][0])<1:
                            if math.floor(tmp[i][0])<tmp[i][0]:
                                if tmp[i] not in tmp1:
                                    tmp1.append(tmp[i])
                            elif math.floor(tmp[i-1][0])<tmp[i-1][0]:
                                if tmp[i-1] not in tmp1:
                                    tmp1.append(tmp[i-1])
                else:
                    print(gt)
                if len(tmp1)>0:
                    tmp2=sorted(tmp1, key=lambda tmp1:tmp1[1], reverse=True)
                    if len(tmp2)>=2:
                        if tmp2[0][0]==temp_gt1:
                            if tmp2[1][0]-tmp2[0][0]>6 and float(tmp2[1][1])/tmp2[0][1]>0.01:
                                temp_gt2=tmp2[1][0]
                                gtstr2=tmp2[1][2]
                            elif tmp2[1][0]-tmp2[0][0]>1 and float(tmp2[1][1])/tmp2[0][1]>0.05:
                                temp_gt2=tmp2[1][0]
                                gtstr2=tmp2[1][2]
                            elif abs(tmp2[1][0]-tmp2[0][0])<1 and float(tmp2[1][1])/tmp2[0][1]>0.001:
                                if len(tmp2)>2 and tmp2[1][1]>tmp2[2][1]:
                                    temp_gt2=tmp2[1][0]
                                    gtstr2=tmp2[1][2]
                                elif len(tmp2)==2:
                                    temp_gt2=tmp2[1][0]
                                    gtstr2=tmp2[1][2]
                        else:
                            print("not eq",gt,temp_gt1,tmp2)
                            print(sorted(site_gts, key=lambda site_gts:site_gts[0], reverse=True))
                            print("all",site_gts)
                    elif len(tmp1)==1:
                        if tmp2[0][0]!=temp_gt1:
                            print("=1",gt,tmp1)
                            print("all",site_gts)
                
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s;%s\n"%(gt,temp_gt1,temp_gt2,",".join(rate),site_gts[0][1],";".join(depthlist),gtstr1,gtstr2))
    else:
        print("other",gt,site_gts)
        
        
infile.close()
outfile.close()
endtime=datetime.datetime.now()
print(time.strftime("%y/%m/%d %H:%M:%S",time.localtime(time.time())))
print("time used :",str(endtime-starttime))
print("well done!")    
