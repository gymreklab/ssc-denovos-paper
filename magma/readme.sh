#!/bin/bash

# Download ref files (b37, since they only have refdata for that)
wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip # Ref data, b37
unzip g1000_eur.zip
wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip # gene locations, b37
unzip NCBI37.3.zip

################ Preprocess: Get TR sets #################
# Get our locations list and convert to b37
./extract_locs.py 2 -1 1 | awk '{print "chr"$0}' > caseonly_mutations_hg38.tab
./extract_locs.py 1 -1 1 | awk '{print "chr"$0}' > ctrlonly_mutations_hg38.tab
./extract_locs.py 2 -1 0.01 | awk '{print "chr"$0}' > caseonly_mutations_hg38_rare.tab
./extract_locs.py 1 -1 0.01 | awk '{print "chr"$0}' > ctrlonly_mutations_hg38_rare.tab
./extract_locs.py 2 0.05 1 | awk '{print "chr"$0}' > caseonly_mutations_hg38_common.tab
./extract_locs.py 1 0.05 1 | awk '{print "chr"$0}' > ctrlonly_mutations_hg38_common.tab

cat /storage/mgymrek/gtex/annotations/coding.bed  /storage/mgymrek/gtex/annotations/5utr.bed | sed 's/chr//' > coding_utr_hg19.bed

for cc in case ctrl
do
    for suffix in "" "_rare" "_common"
    do
	liftOver ${cc}only_mutations_hg38${suffix}.tab /storage/resources/dbase/human/hg19/hg38ToHg19.over.chain ${cc}only_mutations_hg19${suffix}.tab unMapped
	cat ${cc}only_mutations_hg19${suffix}.tab | awk '{print $1"_"$2 "\t" $1 "\t" $2}' | sed 's/chr//g' > ${cc}only_mutations_hg19_snplocs${suffix}.tab

	# All of them
	magma --annotate --snp-loc ${cc}only_mutations_hg19_snplocs${suffix}.tab --gene-loc NCBI37.3.gene.loc --out ${cc}${suffix} > /dev/null
	cat ${cc}${suffix}.genes.annot | grep -v "^#" | awk '{print $1 " "}' | tr -d '\n' | awk -v prefix=${cc}${suffix} '{print prefix " " $0}'

	# Only exon (coding/UTR)
	#cat ${cc}only_mutations_hg19_snplocs${suffix}.tab | \
	#    awk '{print $2 "\t" $3 "\t" $3+1 "\t" $0}' | \
	#    intersectBed -a stdin -b coding_utr_hg19.bed -wa -wb | cut -f 4-6 | sort | uniq > ${cc}only_mutations_hg19_snplocs${suffix}_exon.tab
	#magma --annotate --snp-loc ${cc}only_mutations_hg19_snplocs${suffix}_exon.tab  --gene-loc NCBI37.3.gene.loc --out ${cc}${suffix}_exon > /dev/null
	#cat ${cc}${suffix}_exon.genes.annot | grep -v "^#" | awk '{print $1 " "}' | tr -d '\n' | awk -v prefix=${cc}${suffix}_exon '{print prefix " " $0}'
    done
done > gene_annot.tab


################ Preprocess: GWAS gene analysis #################

EA_GWAS=/storage/ileena/denovos5/EA_GWAS/GWAS_EA_excl23andMe.txt
EA_N=1131881

cat /storage/ileena/denovos5/ASD_GWAS_iPSYC_PGC/iPSYCH-PGC_ASD_Nov2017 | awk '{print $2 "\t" $1 "\t" $3 "\t" $9}' > iPSYCH-PGC_ASD_Nov2017_snplocs.txt
ASD_GWAS=iPSYCH-PGC_ASD_Nov2017_snplocs.txt
ASD_N=46351

echo "SNP,CHR,POS,P" | sed 's/,/\t/g' > ckqny.scz2snpres.snplocs.txt
cat /storage/s1saini/pgc_analysis/ckqny.scz2snpres | grep -v hg19chr | awk '{print $2 "\t" $1 "\t" $5 "\t" $9}' | sed 's/chr//g' >> ckqny.scz2snpres.snplocs.txt
SCZ_GWAS=ckqny.scz2snpres.snplocs.txt
SCZ_N=79845 # TODO: Ileena can you check? I got from the paper 34241+45604

for TRAIT in EA ASD SCZ
do
    ANNOT_PREFIX=${TRAIT}
    GENE_PREFIX=${TRAIT}_GENE
    GWAS_FILE=
    GWAS_N=
    if [ "x$TRAIT" == "xEA" ]; then
	GWAS_FILE=${EA_GWAS}
	GWAS_N=${EA_N}
    fi
    if [ "x$TRAIT" == "xASD" ]; then
	GWAS_FILE=${ASD_GWAS}
	GWAS_N=${ASD_N}
    fi
    if [ "x$TRAIT" == "xSCZ" ]; then
	GWAS_FILE=${SCZ_GWAS}
	GWAS_N=${SCZ_N}
    fi

    # MAGMA step 1 - annotate genes from GWAS
    magma --annotate --snp-loc ${GWAS_FILE} --gene-loc NCBI37.3.gene.loc --out ${ANNOT_PREFIX}
    # MAGMA step 2 - gene analysis
    magma --bfile g1000_eur --gene-annot ${ANNOT_PREFIX}.genes.annot --out ${GENE_PREFIX}  --pval ${GWAS_FILE} N=${GWAS_N}
done

################ Gene set analysis #################
for TRAIT in EA ASD SCZ
do
    magma --gene-results ${TRAIT}_GENE.genes.raw --set-annot gene_annot.tab --out ${TRAIT}_GS
done

# print results to terminal
for trait in ASD SCZ EA; do echo "##### ${trait} #####"; cat ${trait}_GS.gsa.out; done
