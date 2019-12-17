#-----------------------------------
# Authors: a.lobley@qmul.ac.uk [hfx532]
#-----------------------------------
# 
#-----------------------------------
#$ -cwd
#$ -l h_rt=10:00:00,h_vmem=20G
#$ -N applyBFilt
#$ -tc 10
#$ -t 1-499
#-----------------------------------
module load vcftools
module load samtools
module load gatk
module load bedtools
module load star/2.6.1a
#-----------------------------------
# reassign -RMQT all to be 60
# allow u cigar string
# check mark dullicates
#-----------------------------------
KNOWN="/data/BCI-BioInformatics/anna/genomes/ref/1000G_phase1.snps.high_confidence.hg38.vcf"
MILLS="/data/BCI-BioInformatics/anna/genomes/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf"
OMNI="/data/BCI-BioInformatics/anna/datasets/1000G_omni2.5.hg38.vcf.gz"
HAPMAP="/data/BCI-BioInformatics/anna/datasets/hapmap_3.3.hg38.vcf.gz"
DBSNP="/data/BCI-BioInformatics/anna/dbsnp/dbsnp138.vcf"
GATK="/share/apps/centos7/gatk/4.0.8.1/bin/gatk"
GATK_TMP="."
SAMTOOLS=`which samtools`
VCFTOOLS=`which vcftools`
STAR_REF="/data/BCI-BioInformatics/DATASETS/STAR_REF_HG38/Genome"
genome="/data/BCI-BioInformatics/DATASETS/STAR_REF_HG38/hg38.fa"
#----------------------------------------------------------------
#  GATK PROCESS
#----------------------------------------------------------------
#
#  First sort bam files  and extract fastq
#----------------------------------------------------------------
samtools sort -T $SGE_TASK_ID -n $SGE_TASK_ID/$SGE_TASK_ID.bam -o $SGE_TASK_ID/$SGE_TASK_ID.sorted.bam
#----------------------------------------------------------------
if [[ ! -f $SGE_TASK_ID/$SGE_TASK_ID\_r1.fq ]]
then
bedtools bamtofastq  -i $SGE_TASK_ID/$SGE_TASK_ID.sorted.bam \
                     -fq   $SGE_TASK_ID/$SGE_TASK_ID\_r1.fq \
                     -fq2  $SGE_TASK_ID/$SGE_TASK_ID\_r2.fq
fi
date
#----------------------------------------------------------------
# extract header line for re-use
#----------------------------------------------------------------
samtools view -h $SGE_TASK_ID/$SGE_TASK_ID.sorted.bam | grep ^\@ > $SGE_TASK_ID.head
#----------------------------------------------------------------
#rm -rf $SGE_TASK_ID.download/*bam
#rm -rf $SGE_TASK_ID.download/*.bai
ID=`grep \@RG $SGE_TASK_ID.head | cut -f 2 `
LB=`grep \@RG $SGE_TASK_ID.head | cut -f 2 | sed s/ID..//g | cut -d "_" -f 1,2,3`
PL="Illumina"
PU=`grep \@RG $SGE_TASK_ID.head | cut -f 2 | sed s/ID..//g | cut -d "_" -f 4,5,6 `
SM=`grep \@RG $SGE_TASK_ID.head | cut -f 3 | sed s/SM.//g`
#----------------------------------------------------------------
ln --force -s $STAR_REF/* $SGE_TASK_ID/
cd $SGE_TASK_ID
#----------------------------------------------------------------
FQ1=$SGE_TASK_ID"_r1.fq"
FQ2=$SGE_TASK_ID"_r2.fq"
#----------------------------------------------------------------
# Make alignment to genome with v permissive threshold
##----------------------------------------------------------------
STAR --runMode alignReads --readFilesType Fastx \
        --genomeDir . \
        --readFilesCommand cat \
        --readFilesIn $FQ1,$FQ2  \
        --outFileNamePrefix $SGE_TASK_ID \
        --runThreadN 1 \
        --outFilterMultimapScoreRange 1 \
        --outFilterMultimapNmax 999 \
        --outFilterMismatchNmax 20 \
        --alignIntronMax 100000 \
        --alignMatesGapMax 100000 \
        --sjdbScore 2 \
        --alignSJDBoverhangMin 2 \
        --outFilterMatchNmin 25 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMismatchNoverLmax 0.3 \
        --outSJfilterOverhangMin 6 6 6 6 \
        --outSJfilterCountUniqueMin 1 1 1 1 \
        --outSJfilterCountTotalMin 1 1 1 1 \
        --outSAMstrandField intronMotif \
        --outSAMmultNmax 999 \
        --outSAMattributes All \
        --outSAMattrRGline $ID LB:$LB PL:$PL PU:$PU SM:$SM \
        --limitSjdbInsertNsj 900000 \
        --outReadsUnmapped None \
        --twopassMode Basic \
        --outSAMtype BAM  SortedByCoordinate \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimSegmentReadGapMax 3 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimOutJunctionFormat 1
#----------------------------------------------------------------
# mv data to aligned set
##----------------------------------------------------------------
mv Aligned.sortedByCoord.bam ../$SGE_TASK_ID.bam
#----------------------------------------------------------------
if [[ -f $SGE_TASK_ID/$SGE_TASK_ID\Aligned.sortedByCoord.out.bam ]]
then
    mv $SGE_TASK_ID/$SGE_TASK_ID\Aligned.sortedByCoord.out.bam $SGE_TASK_ID.bam
    $SAMTOOLS index $SGE_TASK_ID.bam
fi
#-----------------------------------
if [[ ! -f "$SGE_TASK_ID.head" ]]
then
    samtools view -h $SGE_TASL_ID.bam | grep ^\@RG > $i.head
fi
#-----------------------------------
# Piccard Mark Duplicates
##----------------------------------------------------------------
  $GATK MarkDuplicates -I $SGE_TASK_ID.bam \
                       -M $SGE_TASK_ID.metrics \
                       -O $SGE_TASK_ID.ubam
#-----------------------------------
# SplitNCigarReads and reassign mapping qualities
#-----------------------------------
  $GATK  SplitNCigarReads \
         -R $genome -I $SGE_TASK_ID.ubam \
         -O $SGE_TASK_ID.split.bam \
         -RF MappingQualityReadFilter \
         -RF CigarContainsNoNOperator \
         -RF MappingQualityNotZeroReadFilter \
         --minimum-mapping-quality 2 \
         --refactor-cigar-string=T \
         --maximum-mapping-quality 255 \
         --skip-mapping-quality-transform=F   \
         --TMP_DIR=$GATK_TMP
#-----------------------------------
# Indel Realignment and Base Recalibration
#-----------------------------------
  $GATK  BaseRecalibrator \
          --bqsr-baq-gap-open-penalty 30 \
          --known-sites ${KNOWN} \
          --known-sites ${MILLS} \
          -R ${genome} -I $SGE_TASK_ID.split.bam \
          --TMP_DIR $GATK_TMP \
          -O $SGE_TASK_ID.final.rnaseq.grp
#-----------------------------------
#        --create-output-bam $SGE_TASK_ID.final.bam
#
#($SAMTOOLS view -H $SGE_TASK_ID.final.bam; $SAMTOOLS view $SGE_TASK_ID.final.bam | grep -w 'NH:i:1') \
#|$SAMTOOLS view -Sb -  > $SGE_TASK_ID.final.uniq.bam
#-----------------------------------
  $GATK ApplyBQSR -R ${genome} -O $SGE_TASK_ID.final.bam \
                    --create-output-bam-md5=T \
                    --create-output-bam-index \
                    -I $SGE_TASK_ID.split.bam \
                    -bqsr $SGE_TASK_ID.final.rnaseq.grp
#-----------------------------------
#  ($SAMTOOLS view -H $SGE_TASK_ID.final.bam; $SAMTOOLS view $SGE_TASK_ID.final.bam | grep -w 'NH:i:1') | $SAMTOOLS view -Sb - > $SGE_TASK_ID.uniq.bam
#-----------------------------------
#  rm -f $SGE_TASK_ID.split.bam
#  rm -f $SGE_TASK_ID.gbam
#-----------------------------------
  $SAMTOOLS index $SGE_TASK_ID.split.bam
#-----------------------------------
  # fix absolute path in dict file
 # sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  #echo "${bam.join('\n')}" > $SGE_TASK_ID.bam.list
#-----------------------------------
  # Variant calling
  $GATK  HaplotypeCaller \
          -R $genome -I $SGE_TASK_ID.final.bam \
          --dont-use-soft-clipped-bases \
          --emit-ref-confidence GVCF \
          --standard-min-confidence-threshold-for-calling 10.0 \
          -O $SGE_TASK_ID.gatk4.vcf.gz
#-------------------------------------------
# Steps do not work just yet, syntax problem(s)
#-------------------------------------------
   $GATK VariantRecalibrator \
          -an QD -an MQRankSum -an ReadPosRankSum \
          -an FS -an MQ -mode SNP \
          -an SOR -an InbreedingCoeff \
          -resource:1000G,known=false,training=true,truth=false,prior=10.0,$KNOWN \
          -resource:hapmap,known=false,training=true,truth=true,prior=15.0,$HAPMAP \
          -resource:omni,known=false,training=true,truth=false,prior=12.0,$OMNI \
          -resource:dbsnp,known=true,training=false,truth=false,prior=6.0,$DBSNP \
          --variant input $SGE_TASK_ID.gatk4.vcf.gz \
          -recal-file $SGE_TASK_ID.recal \
          -tranches-file $SGE_TASK_ID.snv.tranches \
          -rscript-file  $SGE_TASK_ID.snv.plots.R \
          -R $genome --max-gaussians 4
#------------------------------------
# Apply recalibration(s) to snps
#------------------------------------
  $GATK  ApplyCalibration -i $SGE_TASK_ID.gatk4.vcf \
                   -V $SGE_TASK_ID.recal \
                   -R $genome \
                   --ts_filter_level 99.0 \
                   -mode SNP \
                   -tranchesFile $SGE_TASK_ID.snv.tranches \
                   -o $SGE_TASK_ID.filt.vcf.gz
#-----------------------------------
# Apply filter to variant snps
#-----------------------------------
  $GATK  VariantFiltration \
         -R $genome \
         -selectMode SNP \
         -V $SGE_TASK_ID.filt.vcf.gz \
         -O $SGE_TASK_ID.final-gsnps.vcf.gz \
         --cluster-size 3 \
         --cluster-window-size 35  \
         --filter-name DP --filter-expression " DP < 2.0 " \
         --filter-name QD --filter-expression " QD < 2.0 " \
         --filter-name FD --filter-expression " FS > 60.0 " \
         --filter-name MQ --filter-expression " MQ < 30.0 " \
         --filter-name MQRankSum  --filter-expression "MQRankSum < -12.5 " \
         --filter-name ReadPosRankSum --filter-expression "ReadPosRankSum < -8.0"
#-----------------------------------
# Apply filter to variant indels
#-----------------------------------
$GATK  VariantFiltration \
         -R $genome \
         -selectMode INDELS \
         -V $SGE_TASK_ID.filt.vcf.gz \
         -O $SGE_TASK_ID.final-indels.vcf.gz \
         --cluster-size 3 \
         --cluster-window-size 35  \
         --filter-name DP --filter-expression " DP < 2.0 " \
         --filter-name QD --filter-expression " QD < 2.0 " \
         --filter-name FD --filter-expression " FS > 60.0 " \
         --filter-name MQ --filter-expression " MQ < 30.0 " \
         --filter-name MQRankSum  --filter-expression "MQRankSum < -12.5 " \
         --filter-name ReadPosRankSum --filter-expression "ReadPosRankSum < -8.0"
#-----------------------------------
# Extract variations to output i.e. convert gvcf to vcf
#-----------------------------------
$GVCF/bin/extract_variants < $SGE_TASK_ID.final-gsnps.vcf.gz | bgzip  > $SGE_TASK_ID.final-snps.vcf
#-----------------------------------