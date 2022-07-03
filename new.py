
import sys
import numpy as np
import gzip


def whole_vcf(thresh):
    vcf=open("am.vcf","r")
    remlist=[]
    het_dp=[]
    total_dp=[]
    hetcount=0
    totalcount=0
    while True:
        line=vcf.readline()
        if line=="":
            break
        if line[0]=="#":
            continue
        totalcount+=1
        line=line.split("\t")
        loc_id=line[2]
        chr=line[0]
        pos=line[1]
        gt=[x.split(":")[0] for x in line[9::]]
        het_count=gt.count("0/1")+gt.count("1/0")
        ref_count=gt.count("0/0")
        alt_count=gt.count("1/1")
        gt_total=len(gt)-gt.count("./.")
        het_ratio=het_count/float(gt_total)
        ref_ratio=ref_count/float(gt_total)
        alt_ratio=alt_count/float(gt_total)
        dp=[x.split(":")[3] for x in line[9::]]
        while "." in dp:
            dp.remove(".")
        for i,d in enumerate(dp):
            dp[i]=int(d)
        dp_site_mean=np.nanmean(dp)
        if het_ratio>=thresh:
            hetcount+=1
            het_dp.append(dp_site_mean)
            remlist.append((chr,pos))
        else:
            total_dp.append(dp_site_mean)
    print("Depth in sites >={0}% heterozygous:".format(thresh*100),np.mean(het_dp))
    print("Depth in other sites:",np.nanmean(total_dp))
    print("% of >{0}% heterozyous sites:".format(thresh*100),hetcount*100/float(totalcount))
    vcf.close()
    return remlist


def write_out(remlist):
    out= open("n.vcf","w")
    vcf= open("am.vcf","r")
    for line in vcf:
        if line[0]=="#":
            out.write(line)
        else:
            chr=line.split("\t")[0]
            pos=line.split("\t")[1]
            if (chr,pos) not in remlist:
                out.write(line)

remlist=whole_vcf(0.2)
remlist=set(remlist)
write_out(remlist)