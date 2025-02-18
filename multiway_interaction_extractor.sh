
## filter for UU interactions on chr2 with mapq > 1 from pairs file
infile=$1

prefix=$2 #runid
name=`basename $infile .pairs.gz`

probeBed=$3

outfolder=$4

pairtools select '(pair_type=="UU") and (chrom1=="chr2") and (chrom2=="chr2") and (mapq1>1) and (mapq2 > 1)' -o ${prefix}_${name}_chr2.pairs.gz ${infile}

## format the pairs file to simple version bedpe file. 10 columns: chrom1 start1 end1 chrom2 start2 end2 strand1 strand2 placeholder readid.
zcat ${prefix}_${name}_chr2.pairs.gz |grep -v "#"| awk '{OFS="\t"; print $1,$2,$13,$15,$4,$14,$16,$6,$7}'| awk '{OFS="\t";if ($8=="-" ){$10=$3;$3=$4;$4=$10}if($9=="-"){$10=$6;$6=$7;$7=$10} print$2,$3,$4,$5,$6,$7,$8,$9,"0",$1}' > ${outfolder}/${prefix}_${name}_chr2.simple.pairs.bedpe

## filter reads with both fragments overlap with probe regions. 
pairToBed -a ${outfolder}/${prefix}_${name}_chr2.simple.pairs.bedpe -b ${probeBed} -type both > ${outfolder}/${prefix}_${name}_chr2_probe_overlap_bothfrag.txt 

## filter reads with fragments overlap with different probe regions
python multiway.py ${outfolder}/${prefix}_${name}_chr2_probe_overlap_bothfrag.txt ${outfolder}/${prefix}_${name}_chr2_probe_overlap_bothfrag ${outfolder}/${prefix}

###usage: sh multiway_interaction_extractor.sh Run15_Hudep2_E1KO34_NlaIII_pro_CTCF/pairs/*pairs.gz Run15 run3_6_peaks2.bed 20240417_result 
