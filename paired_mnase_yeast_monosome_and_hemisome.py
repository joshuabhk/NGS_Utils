from pysam import Samfile
import numpy as np

def get_fragment_counts( bam, chrom="chr1", rfac=1, smooth=True ):
    size = bam.lengths[bam.references.index(chrom)]
    
    monosomes = np.zeros(ceil(size/rfac))
    hemisomes = np.zeros(ceil(size/rfac))
    unmatched = {}
    
    #count monosomes and hemisomes for a debugging purpose
    mcount, hcount = 0, 0 

    mapq_cutoff = 20
    for aln in bam.fetch(reference=chrom):
        if aln.is_supplementary :
            continue

        if aln.query_name not in unmatched:
            unmatched[aln.query_name] = aln
            continue

        aln2 = unmatched.pop(aln.query_name)

        if aln.mapq < mapq_cutoff and aln2.mapq < mapq_cutoff :
            continue

        assert aln.reference_name == aln2.reference_name

        positions = [aln.reference_start, aln.reference_end, 
                     aln2.reference_start, aln2.reference_end ]
        l1, l2 = min(positions), max(positions)
        insert = l2-l1

        #focal point
        loc = (l1 + l2) / 2
        rloc = int(loc/rfac)
        
        rbegin = int(l1/rfac)
        rend = ceil(l2/rfac)

        if 130 < insert < 170 :
            mcount += 1
            #print(insert, rloc)
            if smooth:
                delta = 2/(rend-rbegin)
                monosomes[rbegin:rloc] += [(i+1)*delta for i in range(rloc-rbegin)]
                monosomes[rloc:rend] += [(i+1)*delta for i in reversed(range(rend-rloc))]

            else: 
                monosomes[rbegin:rend] += 1
            
        if 55 < insert < 95 :
            hcount += 1
            if smooth :
                delta = 2/(rend-rbegin)
                hemisomes[rbegin:rloc] += [(i+1)*delta for i in range(rloc-rbegin)]
                hemisomes[rloc:rend] += [(i+1)*delta for i in reversed(range(rend-rloc))]

            else:
                hemisomes[rbegin:rend] += 1
            
            
    print( chrom,":", mcount, hcount )        

    #total = sum(monosomes)+sum(hemisomes)

    #monosomes = monosomes*1000000/monosomes.sum()
    #hemisomes = hemisomes*1000000/hemisomes.sum()
    
    return monosomes, hemisomes


def write_monosome_hemisome_bedgraphs(bamfn, 
                                      monosome_label="mono",
                                      hemisome_label="hemi",
                                      postfix="bedgraph",
                                      smooth=False
                                     ) :
    bam = Samfile(bamfn)
    print(bam.lengths)
    print(bam.references)
    outbase = bamfn.replace( ".bam", "" )
    
    if smooth == True :
        postfix = "smooth."+postfix
        
    mfn = ".".join([outbase, monosome_label, postfix])
    hfn = ".".join([outbase, hemisome_label, postfix])
    print(mfn, hfn)
    
    mfp = open( mfn, 'w')
    hfp = open( hfn, "w")
    
    for chrom in bam.references :
        monosomes, hemisomes = get_fragment_counts( bam, 
                                                   chrom=chrom, 
                                                   smooth=smooth )
        
        for i, (j,k) in enumerate(zip(monosomes,hemisomes)) :
            print(chrom, i, i+1, j, file=mfp, sep="\t")
            print(chrom, i, i+1, k, file=hfp, sep="\t")

import sys
if __name__ == "__main__" :
    bamfns = sys.argv[1:]
    for bamfn in bamfns:
        write_monosome_hemisome_bedgraphs(bamfn)
