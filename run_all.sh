# #!/bin/bash

# ############# Part I: Setup and Initialization
# # Paths and directories
# SAMPLES=/home/mkonczal/Teaching/GEiP/Data/SamplesLibraries.txt
# FASTQ_DIR=/home/mkonczal/Teaching/GEiP/Data/Fastq
# REF=/home/mkonczal/Teaching/GEiP/Data/Reference/SBS_final.scaffolds.fasta

# # Create directories
# mkdir -p GEiP/Lab1
# cd GEiP/Lab1

# # Extract scaffold1 from the reference genome
# seqkit grep -p "scaffold1" $REF > scaffold1.fasta
# samtools faidx scaffold1.fasta
# cut -f 2 scaffold1.fasta.fai

# # Function to process each individual
# process_individual() {
#     local INDIVIDUAL=$1
#     local RUN1=$2
#     local RUN2=$3

#     echo "Processing individual: $INDIVIDUAL"

#     # Copy FASTQ files
#     cp ${FASTQ_DIR}/${RUN1}/${RUN1}_pass_1.fastq.gz .
#     cp ${FASTQ_DIR}/${RUN1}/${RUN1}_pass_2.fastq.gz .
#     cp ${FASTQ_DIR}/${RUN2}/${RUN2}_pass_1.fastq.gz .
#     cp ${FASTQ_DIR}/${RUN2}/${RUN2}_pass_2.fastq.gz .

#     # Define file variables
#     f1=${RUN1}_pass_1.fastq.gz
#     f2=${RUN1}_pass_2.fastq.gz
#     f3=${RUN2}_pass_1.fastq.gz
#     f4=${RUN2}_pass_2.fastq.gz

#     # Quality control and filtering
#     fastp -i $f1 -I $f2 -o ${RUN1}_pass_1.filt.fastq -O ${RUN1}_pass_2.filt.fastq
#     fastp -i $f3 -I $f4 -o ${RUN2}_pass_1.filt.fastq -O ${RUN2}_pass_2.filt.fastq

#     fastqc ${RUN1}_pass_1.filt.fastq ${RUN1}_pass_2.filt.fastq
#     fastqc ${RUN2}_pass_1.filt.fastq ${RUN2}_pass_2.filt.fastq

#     # Align reads to the reference genome
#     bwa index scaffold1.fasta

#     bwa mem -t 10 -R "@RG\tID:${RUN1}\tSM:${INDIVIDUAL}\tLB:${RUN1}\tPL:ILLUMINA\tPU:lib1_unit" scaffold1.fasta ${RUN1}_pass_1.filt.fastq ${RUN1}_pass_2.filt.fastq | samtools view -F 4 -o ${RUN1}.Mapped.bam
#     bwa mem -t 10 -R "@RG\tID:${RUN2}\tSM:${INDIVIDUAL}\tLB:${RUN2}\tPL:ILLUMINA\tPU:lib1_unit" scaffold1.fasta ${RUN2}_pass_1.filt.fastq ${RUN2}_pass_2.filt.fastq | samtools view -F 4 -o ${RUN2}.Mapped.bam

#     # Clean up intermediate files
#     rm $f1 $f2 $f3 $f4

#     # Sort BAM files
#     samtools sort -T bam ${RUN1}.Mapped.bam > ${RUN1}.Mapped.sorted.bam
#     samtools sort -T bam ${RUN2}.Mapped.bam > ${RUN2}.Mapped.sorted.bam

#     rm ${RUN1}.Mapped.bam ${RUN2}.Mapped.bam

#     # Mark duplicates
#     picard MarkDuplicates REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true AS=true I=${RUN1}.Mapped.sorted.bam M=test.metric_${RUN1}.txt O=${RUN1}.Mapped.sorted_DupRmv.bam 2> MarkDup_${RUN1}.log
#     picard MarkDuplicates REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true AS=true I=${RUN2}.Mapped.sorted.bam M=test.metric_${RUN2}.txt O=${RUN2}.Mapped.sorted_DupRmv.bam 2> MarkDup_${RUN2}.log

#     rm ${RUN1}.Mapped.sorted.bam ${RUN2}.Mapped.sorted.bam

#     # Index BAM files
#     samtools index ${RUN1}.Mapped.sorted_DupRmv.bam
#     samtools index ${RUN2}.Mapped.sorted_DupRmv.bam

#     # Generate flagstats
#     samtools flagstats ${RUN1}.Mapped.sorted_DupRmv.bam
#     samtools flagstats ${RUN2}.Mapped.sorted_DupRmv.bam

#     # Merge BAM files for the individual
#     samtools merge ${INDIVIDUAL}.bam ${RUN1}.Mapped.sorted_DupRmv.bam ${RUN2}.Mapped.sorted_DupRmv.bam
#     samtools index ${INDIVIDUAL}.bam

#     # View alignment
#     # samtools tview --reference scaffold1.fasta ${INDIVIDUAL}.bam
# }

# # Process individuals
# process_individual C_ruf_09 SRR7054135 SRR7054162
# process_individual C_pyg_26 SRR7054133 SRR7054147

# ############# Part II: UCE Analysis
# mkdir -p GEiP/Lab2/scaffold1
# cd /home/st2/GEiP/Lab2

# # Convert fasta to 2bit format
# faToTwoBit ../Lab1/scaffold1.fasta scaffold1.2bit
# twoBitInfo scaffold1.2bit sizes.tab

# ln -s /home/mkonczal/Teaching/GEiP/Data/UCE-probes/uce-5k-probes.fasta .
# samtools faidx uce-5k-probes.fasta
# bwa index uce-5k-probes.fasta

# UCEprobe=/home/mkonczal/Lab2/uce-5k-probes.fasta

# phyluce_probe_run_multiple_lastzs_sqlite --db scaffold1.sqlite --output scaffold1-genome-lastz --scaffoldlist scaffold1 --genome-base-path ./ --probefile ${UCEprobe} --cores 10

# Create scaffold1.conf file
# cat <<EOL > scaffold1.conf
# [scaffolds]
# scaffold1:/home/st2/GEiP/Lab2/scaffold1/scaffold1.2bit
# EOL

# phyluce_probe_slice_sequence_from_genomes --lastz scaffold1-genome-lastz --flank 1000 --output OUT --conf scaffold1.conf --name-pattern uce-5k-probes.fasta_v_scaffold1.lastz.clean

# mkdir UCE_regions

# grep "uce" OUT/scaffold1.fasta | cut -d '|' -f 2,3,4,6 | sed -e 's/|/\t/g' | sed -e 's/contig://g' | sed -e 's/slice://g'| sed -e 's/uce://g' | sed -e 's/orient://g' | sed -e 's/uce-/uce_/g' | sed -e s/"'"//g | sed -e 's/{+}/forward/g' | sed -e 's/{-}/reverse/g'| sed -e 's/-/\t/g' > UCE_regions/scaffold1_UCE_regions.txt 

# mkdir UCE_regions/forward
# mkdir UCE_regions/reverse

# grep 'forward' UCE_regions/scaffold1_UCE_regions.txt | cut -f 1,2,3 > UCE_regions/forward/scaffold1_UCE_forward_orient_regions.txt
# # did not yield anything ^ *shrug*
# grep 'reverse' UCE_regions/scaffold1_UCE_regions.txt | cut -f 1,2,3 > UCE_regions/reverse/scaffold1_UCE_reverse_orient_regions.txt

# UCE=UCE_regions/reverse/scaffold1_UCE_reverse_orient_regions.txt
# REF=../Lab1/scaffold1.fasta
# BAM=../Lab1/C_pyg_26.bam

# samtools mpileup -l ${UCE} -f ${REF} ${BAM} > out.reverse.mpileup
# bcftools mpileup --threads 10 -Ou -Q 30 -q 30 -C 50 -a AD,DP -R ${UCE} -f ${REF} ${BAM} | bcftools call --threads 10 -c -Ob > out.reverse.bcf
# cat <<EOL > bam_list.txt
# /home/st2/GEiP/Lab1/C_pyg_26.bam
# /home/st2/GEiP/Lab1/C_ruf_09.bam
# EOL

# bcftools mpileup --threads 10 -Ou -Q 30 -q 30 -C 50 -a AD,DP -R ${UCE} -f ${REF} -b bam_list.txt | bcftools call --threads 10 -c -Ob > out.reverse_2samples.bcf

# BCF=out.reverse.bcf
# bcftools query -f '%QUAL\t%MQ\t%DP\n' ${BCF} > Stats_QualMQDP.txt
# bcftools stats out.reverse.bcf
# bcftools view -v snps -e 'QUAL < 20 || MQ < 40 || FORMAT/DP < 3 || FORMAT/DP > 100' out.reverse.bcf > out.reverse.filtered.vcf

# UCE=UCE_regions/reverse/scaffold1_UCE_reverse_orient_regions.txt
# REF=../Lab1/scaffold1.fasta
# BAM=../Lab1/C_ruf_09.bam

# samtools mpileup -l ${UCE} -f ${REF} ${BAM} > out.reverse_ruf_09.mpileup
# bcftools mpileup --threads 10 -Ou -Q 30 -q 30 -C 50 -a AD,DP -R ${UCE} -f ${REF} ${BAM} | bcftools call --threads 10 -c -Ob > out.reverse_ruf_09.bcf
# BCF=out.reverse_ruf_09.bcf
# bcftools query -f '%QUAL\t%MQ\t%DP\n' ${BCF} > Stats_QualMQDP_ruf_09.txt
# bcftools stats out.reverse_ruf_09.bcf
# bcftools view -v snps -e 'QUAL < 25 || MQ < 30 || FORMAT/DP < 2 || FORMAT/DP > 100' out.reverse_ruf_09.bcf > out.reverse_ruf_09.filtered.vcf

# ############# Part III: Lab 3
# mkdir ~/GEiP/Lab3 && cd ~/GEiP/Lab3
# ln -s ../Lab1/scaffold1.fasta* .

# galgal=/home/mkonczal/Teaching/GEiP/Data/galGal6
# scaffold=scaffold1.fasta

# blastn -query ${scaffold} -db ${galgal}/Gallus_gallus.GRCg6a.dna_rm.toplevel.fa -outfmt 6 > Scaffold1Chicken.blastout
# cut -f2 Scaffold1Chicken.blastout | sort | uniq -c

# scores=/home/mkonczal/Teaching/GEiP/utilities/HoxD55
# hom_chicken_chr=5

# lastz ${galgal}/split/${hom_chicken_chr}.fa ${scaffold} --ambiguous=iupac --hspthresh=2200 --inner=2000 --ydrop=3400 --gappedthresh=10000 --scores=${scores} --chain --format=axt > bGalGal6_chr${hom_chicken_chr}.axt

# alignment=bGalGal6_chr${hom_chicken_chr}.axt
# chicken_2bit=${galgal}/Gallus_gallus.GRCg6a.dna_rm.toplevel.2bit
# faToTwoBit scaffold1.fasta scaffold1.2bi
# biegus_2bit=scaffold1.2bit
# axtChain -minscore=5000 -linearGap=loose $alignment $chicken_2bit $biegus_2bit bgalgalChr${hom_chicken_chr}_scaff1.chain
# chainSort bgalgalChr${hom_chicken_chr}_scaff1.chain sorted_bgalgalChr${hom_chicken_chr}_scaff1.chain

# grep "chain" sorted_bgalgalChr${hom_chicken_chr}_scaff1.chain | wc -l

# awk '{print "scaffold1\t" $2}' scaffold1.fasta.fai > scaffold1.chrom.size
# ln -s ${galgal}/split/${hom_chicken_chr}.fa
# samtools faidx ${hom_chicken_chr}.fa
# awk -v chr="${hom_chicken_chr}" '{print chr "\t" $2}' ${hom_chicken_chr}.fa.fai > ${hom_chicken_chr}.chrom.size

# chainNet sorted_bgalgalChr${hom_chicken_chr}_scaff1.chain 5.chrom.size scaffold1.chrom.size all.net /dev/null
# netChainSubset all.net sorted_bgalgalChr${hom_chicken_chr}_scaff1.chain galGalChr_${hom_chicken_chr}ToSBS_Scaff1.over.chain

# gzip galGalChr_${hom_chicken_chr}ToSBS_Scaff1.over.chain

# chCADD_dir=/home/mkonczal/Teaching/GEiP/Data/chCADD-scores
# mkdir chCADD && cd chCADD

# cp ${chCADD_dir}/Header.txt .
# cp ${chCADD_dir}/5.txt.gz .
# cat Header.txt > chr5_chCADD.tsv
# zcat 5.txt.gz >> chr5_chCADD.tsv

# awk '{print $1,$2-1,$2,$3,$4,$5}' chr5_chCADD.tsv > chr5_chCADD.1based.bed
# rm chr5_chCADD.tsv

# conda activate crossmap
scaffold=scaffold1
# CrossMap bed galGalChr_${hom_chicken_chr}ToSBS_Scaff1.over.chain.gz chCADD/chr${hom_chicken_chr}_chCADD.1based.bed | grep $scaffold | grep -v "Unmap" | cut -f 3,4,5,6,7,8 > chr${hom_chicken_chr}-SBS_CADD.bed

############# Part IV: Lab 4
# cd ~/GEiP/Lab3
# mkdir ~/GEiP/Lab4 && cd ~/GEiP/Lab4
# mkdir C_pyg_26 C_ruf_09

# Function to process VCF files
process_vcf() {
    local INDIVIDUAL=$1
    local VCF_R=$2

    echo "Processing VCF for individual: $INDIVIDUAL"

    # conda activate crossmap
    cd ${INDIVIDUAL}

    vcf2bed --max-mem 4G < ${VCF_R} > vcf_${INDIVIDUAL}_reverse.bed
    grep -v "INDEL" vcf_${INDIVIDUAL}_reverse.bed > vcf_${INDIVIDUAL}_reverse_indelRm.bed
    CADD=../../Lab3/chr5-SBS_CADD.bed
    bed=vcf_${INDIVIDUAL}_reverse_indelRm.bed
    bedtools intersect -a $CADD -b $bed -wa -wb > vcf_${INDIVIDUAL}_reverse_indelRm_intersect.bed

    script_path=/home/mkonczal/Teaching/GEiP/scripts
    awk -v b=6 -v e=100 -f ${script_path}/Orient_conversion_reverse_mod.awk vcf_${INDIVIDUAL}_reverse_indelRm_intersect.bed > vcf_${INDIVIDUAL}_reverse_indelRm_intersect_OrientedReverse.bed
    awk -v b=6 -v e=100 -f ${script_path}/SNP_check_forward.awk vcf_${INDIVIDUAL}_reverse_indelRm_intersect_OrientedReverse.bed > vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated.bed
    cut -f 23 vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated.bed | sort | uniq -c

    awk -F'\t' '$23 == "SNP_is_ALT_pp=ref"' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated.bed > vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed

    awk -e ' $20 ~ /^0\/0/ ' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed | cut -f 6 | paste -sd+ - | bc

    awk -e ' $20 ~ /^0\/1/ ' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed | cut -f 6 | paste -sd+ - | bc

    awk -e ' $20 ~ /^1\/1/ ' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed | cut -f 6 | paste -sd+ - | bc

    # obciązenie całkowite
    awk -e '$20 ~ /^0\/1/ {print $6 * 0.5}' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed | paste -sd+ - | bc
    # obciążenie zrealizowane
    awk -e '$20 ~ /^0\/1/ {print $6 * 0.1}' vcf_${INDIVIDUAL}_reverse_indelRm_intersect_annotated_SNP_is_alt.bed | paste -sd+ - | bc

    cd ..
}

# Process VCF files for C_pyg_26
VCF1_r=../../Lab2/out.reverse_C_pyg_26.filtered.vcf
process_vcf C_pyg_26 ${VCF1_r}

# # Process VCF files for C_ruf_09
# VCF2_r=../../Lab2/out.reverse_C_ruf_09.filtered.vcf
# process_vcf C_ruf_09 ${VCF2_r}