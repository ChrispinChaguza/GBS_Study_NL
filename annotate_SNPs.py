#!/usr/bin/env python

import os
import sys
from Bio import SeqIO
from progress.bar import FillingSquaresBar
from datetime import datetime

def main():

    snp_pos_file = sys.argv[1]
    annotfiles=[sys.argv[2]]
    out_file_name = sys.argv[3]

    db={}

    for h in annotfiles:
        db[str(h)]=SeqIO.read(h,"genbank")

    blast=[]

    with open(snp_pos_file,"r") as fhandle:
        for i in fhandle:
            blast.append(str(i).strip())

    with open(out_file_name,"w") as fhandle:
        fhandle.write("Position\tGeneStart\tGeneEnd\tGeneName\tPositionInGene\tLocusTag\tStrand\tPositionInGenome\tProduct\n")

        max_prog_bar = len(blast)
        status = FillingSquaresBar('['+str(datetime.now().strftime("%H:%M:%S"))+']', max=max_prog_bar)

        for q,j in enumerate(blast):
            if str(j)=="":
                continue

            for i in annotfiles:
                annot=SeqIO.read(str(i),"genbank")

                if True:
                    check=False
                    outp=""

                    for k in annot.features:
                        if k.type=="CDS":
                            if k.strand==1:
                                start=k.location.start.position
                                end=k.location.end.position
                            else:
                                start=k.location.end.position
                                end=k.location.start.position

                            product=""
                            gene=""
                            locus=""
                            genePos=""
                            strand=k.strand

                            if int(j)>start:
                                genePos=int(j)-start+1
                            else:
                                genePos=start-int(j)+1

                            if 'product' in k.qualifiers:
                                product=k.qualifiers['product'][0]
                            if 'gene' in k.qualifiers:
                                gene=k.qualifiers['gene'][0]
                            if 'locus_tag' in k.qualifiers:
                                locus=k.qualifiers['locus_tag'][0]

                            if ((int(j)>=start) and (int(j)<=end)) or ((int(j)>=end) and (int(j)<=start)):
                                if 'product' in k.qualifiers:
                                    product=k.qualifiers['product'][0]
                                if 'gene' in k.qualifiers:
                                    gene=k.qualifiers['gene'][0]
                                if 'locus_tag' in k.qualifiers:
                                    locus=k.qualifiers['locus_tag'][0]

                                outp=str(j)+"\t"+str(start)+"\t"+str(end)+"\t"+str(gene)+"\t"+str(genePos)+"\t"+str(locus)+"\t"+str(strand)+"\tCoding/genic\t"+str(product)
                                check=True 
                                break
                            else:
                                pass
                    if check:
                        fhandle.write(str(outp)+"\n")
                    else:
                         fhandle.write(str(j)+"\tN/A\tN/A\tN/A\tN/A\tN/A\t"+str(strand)+"\tIntergenic region\tIntergenic region\n")
            status.next()
    status.finish()

if __name__=="__main__":
    main()
