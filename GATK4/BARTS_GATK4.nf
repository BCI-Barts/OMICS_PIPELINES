/*
 * Copyright (c) 2019, BCI-BARTS.
 *
 *   This file is part of 'GATK4': 
 *   A Nextflow pipeline for Variant Calling with RNASeq NGS data
 *
 *   BARTS-NGS is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   BARTS-NGS is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with grandRNA-NF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
  
/* 
 * 'BARTS-NGS' - A Nextflow pipeline for variant calling with NGS data
 * 
 * This pipeline that reproduces steps from the GATK best practics of SNP 
 * calling with RNAseq data procedure from STAR :
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 * 
 * 
 * Anna Lobley
 * Pauline Marie Fourgoux
 */


/*
 * Define the default parameters
 */ 
params.star       = "/usr/bin/STAR"
params.genome     = "$baseDir/data/genome/hg38.fa"
params.variants   = "$baseDir/data/snps/hg38_known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/snps/hg38_blacklist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"
params.gatk       = '/usr/bin/gatk'
params.gatk_tmp   = '/tmp'
params.picard     = `which picard-tools`
params.paired     = 'True'
params.gtf        = '$baseDir/data/gtf/hg38.gtf'
params.readlen    = 100
params.mappability = '$baseDir/data/mappability/$readlen/hg38.bw'
params.grep       = `which grep`
params.samtools   = `which samtools`
params.tabix      = `which tabix`

log.info """\
BARTS-NGS v 1.1
================================
genome   : $params.genome
star     : $params.star
grep     : $params.grep
picard   : $params.picard
gtf      : $params.gtf
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
gatk_tmp : $params.gatk_tmp
reads_paired: $params.paired
knownIso : $params.knownIso
mappability : $params.mappability
readlen  : $params.readlen
tabix    : $params.tabix
vcftools : $params.vcftools
"""

/*
 *  Parse the input parameters
 */
STAR            = params.star
PICARD          = params.picard
BOWTIE          = params.bowtie
GATK            = params.gatk_launch
GATK_TMP        = file(params.gatk_tmp)
genome_file     = file(params.genome)
variants_file   = file(params.variants)
blacklist_file  = file(params.blacklist)
reads_ch        = Channel.fromFilePairs(params.reads)
genome_bigWig   = file(params.mappability)
knownIso        = file(params.knownIso)
gtf             = file(params.gtf)
paired          = params.paired
readlen         = params.readlen
GREP            = params.grep
VCFTOOLS        = params.vcftools
SAMTOOLS        = params.samtools
TABIX	    	    = params.tabix
/**********
 * PART 1: Data preparation
 *
 * Process 1A: Assume fasta genome index for star is already made
 */

process 'extractFastQ' { 
  tag "$genome.baseName"
  
  input: 
      file fastq from genome_file 
 
  output: 
      file "${genome}.fai" into genome_index_ch  
  
  script:
  """
  $SAMTOOLS faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "${genome.baseName}.dict" into genome_dict_ch

  script:
  """
  \$PICARD CreateSequenceDictionary R=$genome O=${genome.baseName}.dict
  """
}


/*
 * Process 1C: Create STAR index file.
 */

process '1C_prepare_star_genome_index' {
  tag "$genome.baseName"

  input:
      file genome from genome_file
  output:
      file "genome_dir" into genome_dir_ch

  script:
  """
  #!/usr/bin/env bash
  
  mkdir genome_dir
  $STAR --runMode genomeGenerate \
        --genomeDir genome_dir \
        --genomeFastaFiles ${genome} \
        --runThreadN ${task.cpus}
  
  """
}


/*
 * Process 1D: Align to genome with STAR 
 */

process '1D_star_genome_align'{

  tag "$genome.baseName"
  
  input:
    file  genome from genom_dir_ch

  output:
    file  "Aligned.sortedByCoordinate.bam" into galign
    
    script:
    """
    #! /usr/bin/env bash
    
    $STAR --genomeFastaFiles ${genome} \
          --runMode \
          --runThreadN $task.cpus \
          --genomeDir  . \
          --outputDir .
    
    """
}

process '1E_whippet_run'{

  tag "$whippet"
  
  input:
    file windex from whippet_index
    file reads from reads_ch
    file genome_fasta from $genome
    file gtf from $gff
    
  output:
    file "psi.gz" into whippet_psi
    file "tpm.gz" into whippet_tpm
    file "isoform.gz" into whippet_iso
    
  script:
  """
  #!/usr/bin/env bash
  
  # vanilla whippet on known gtf
  
  $julia  $whippet/whippet-quant.jl \
      -x  $whippet_index \
      -fq $fq1,$fq2 \
      -o  $whippet_out
  # novel whippet finding new splice forms
  
  $julia  $whippet/whippet-index.jl \
        --gtf   $gtf \
        --fasta $genome_fasta \
        --bam   $gbam  \
        -x      $bindex \
        --suppress-low-tsl
  $julia $whippet/whippet-quant.jl \
          $fq1 \
          $fq2 \
          -x $bindex \
          --biascorrect
  """

}

process '1F_prepare_vcf_file' {
  tag "$variantsFile.baseName"

  input: 
      file variantsFile from variants_file
      file blacklisted from blacklist_file

  output:
      set file("${variantsFile.baseName}.filtered.recode.vcf.gz"), file("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi") into prepared_vcf_ch
  
  script:  
  """
  \$VCFTOOLS --gzvcf $variantsFile -c \
           --exclude-bed ${blacklisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz
  \$TABIX ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}

/*
 *  END OF PART 1
 *********/


rnaseq_mapping_star' {
  tag "$replicateId"
  input: 
      file genome from genome_file 
      file genomeDir from genome_dir_ch
      set replicateId, file(reads) from reads_ch 
  output: 
      set replicateId, file('Aligned.sortedByCoord.out.bam'), file('Aligned.sortedByCoord.out.bam.bai') into aligned_bam_ch
  script:    
  """
  # ngs-nf-dev Align reads to genome
  # Final read alignments
  #--------------------------------------
  \$STAR --genomeDir genomeDir \
       --readFilesIn $reads \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --two-pass-mode \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:$replicateId LB:library PL:illumina PU:machine SM:GM12878
/*
 *  END OF PART 2
 ******/
/**********
 * PART 3: GATK Prepare Mapped Reads
 *
 * Process 3: Split reads that contain Ns in their CIGAR string.
 *            Creates k+1 new reads (where k is the number of N cigar elements) 
 *            that correspond to the segments of the original read beside/between 
 *            the splicing events represented by the Ns in the original CIGAR.
 */
process '3_rnaseq_gatk_splitNcigar' {
  tag "$replicateId"
  
  input: 
      file genome from genome_file 
      file index from genome_index_ch
      file genome_dict from genome_dict_ch
      set replicateId, file(bam), file(index) from aligned_bam_ch
  output:
      set replicateId, file('split.bam'), file('split.bai') into splitted_bam_ch
  
  script:
  """
  # SplitNCigarReads and reassign mapping qualities
  $GATK  SplitNCigarReads \
          -R $genome -I $bam \
          -O split.bam \
          -RF MappingQualityReadFilter 
          -RF CigarContainsNoNOperator
          -RF MappinqQualityNotZeroReadFilter
          --minimum-mapping-quality 60 
          --refactor-cigar-string=T
          --maximum-mapping-quality 255
          --skip-mapping-quality-transform=F
          --TMP_DIR=$GATK_TMP
  """
          
}
/*
 *  END OF PART 3
 ******/
/***********
 * PART 4: GATK Base Quality Score Recalibration Workflow
 *
 * Process 4: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */
process '4_rnaseq_gatk_recalibrate' {
  tag "$replicateId"
    
  input: 
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set replicateId, file(bam), file(index) from splitted_bam_ch
      set file(variants_file), file(variants_file_index) from prepared_vcf_ch
  output:
      set sampleId, file("${replicateId}.final.uniq.bam"), file("${replicateId}.final.uniq.bam.bai") into (final_output_ch, bam_for_ASE_ch)
  
  script: 
  sampleId = replicateId.replaceAll(/[12]$/,'')
  """
  # Indel Realignment and Base Recalibration
  $GATK  BaseRecalibrator \
          --bsqr-baq-gap-open-penalty=30 \
          -knownSites ${variants_file} \
          -R ${genome} -I ${bam} \
          --TMP_DIR $GATK_TMP \
          -O final.rnaseq.grp \
          --create-output-bam final.bam \
          --create-output-bam-md5=T \
          --create-output-bam-index
  # Select only unique alignments, no multimaps
  ($SAMTOOLS view -H final.bam; $SAMTOOLS view final.bam| grep -w 'NH:i:1') \
  |$SAMTOOLS view -Sb -  > ${replicateId}.final.uniq.bam
  # Index BAM files
  $SAMTOOLS index ${replicateId}.final.uniq.bam
  """
}
/*
 *  END OF PART 4
 ******/
/***********
 * PART 5: GATK Variant Calling
 *
 * Process 5: Call variants with GATK HaplotypeCaller.
 *            Calls SNPs and indels simultaneously via local de-novo assembly of 
 *            haplotypes in an active region.
 *            Filter called variants with GATK VariantFiltration.    
 */
process '5_rnaseq_call_variants' {
  tag "$sampleId"
  input:
      file genome from genome_file
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(bam), file(bai) from final_output_ch.groupTuple()
  
  output: 
      set sampleId, file('final.vcf') into vcf_files
  script:
  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  echo "${bam.join('\n')}" > bam.list
  
  # Variant calling
  $GATK  HaplotypeCaller \
          -R $genome -I bam.list \
          --dont-use-soft-clipped-bases \
          --emit-ref-confidence GVCF \
          --standard-min-confidence-threshold-for-calling 20.0 \
          -O output.gatk.vcf.gz 

  $GATK  VariantFiltration \
         -R $genome \
         -V output.gatk.vcf.gz \
         -O final.vcf \
         --cluster-size 3 \
         --cluster-window-size 35  \
         --filterName FS -filter "FS > 30.0" \
         --filterName QD -filter "QD < 2.0" 
  """
}
/*
 *  END OF PART 5
 ******/
  
/***********
 * PART 6: Post-process variants file and prepare for Allele-Specific Expression and RNA Editing Analysis
 *
 * Process 6A: Post-process the VCF result  
 */
process '6A_post_process_vcf' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  
  input:
      set sampleId, file('final.vcf') from vcf_files
      set file('filtered.recode.vcf.gz'), file('filtered.recode.vcf.gz.tbi') from prepared_vcf_ch 
  output: 
      set sampleId, file('final.vcf'), file('commonSNPs.diff.sites_in_files') into vcf_and_snps_ch
  
  script:
  '''
  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf
  
  $VCFTOOLS --vcf result.DP8.vcf --gzdiff filtered.recode.vcf.gz  --diff-site --out commonSNPs
  '''
}
/* 
 * Process 6B: Prepare variants file for allele specific expression (ASE) analysis
 */
process '6B_prepare_vcf_for_ase' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  
  input: 
      set sampleId, file('final.vcf'), file('commonSNPs.diff.sites_in_files') from vcf_and_snps_ch
  output: 
      set sampleId, file('known_snps.vcf') into vcf_for_ASE
      file('AF.histogram.pdf') into gghist_pdfs
  script:
  '''
  awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed
    
  $VCFTOOLS --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf
  $GREP -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
               |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
               @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
               {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
               >AF.4R
  gghist.R -i AF.4R -o AF.histogram.pdf
  '''
}
/* 
 * Group data for allele-specific expression.
 * 
 * The `bam_for_ASE_ch` emites tuples having the following structure, holding the final BAM/BAI files:
 *  
 *   ( sample_id, file_bam, file_bai )
 * 
 * The `vcf_for_ASE` channel emits tuples having the following structure, holding the VCF file:
 *  
 *   ( sample_id, output.vcf ) 
 * 
 * The BAMs are grouped together and merged with VCFs having the same sample id. Finally 
 * it creates a channel named `grouped_vcf_bam_bai_ch` emitting the following tuples: 
 *  
 *   ( sample_id, file_vcf, List[file_bam], List[file_bai] )
 */
bam_for_ASE_ch
    .groupTuple()
    .phase(vcf_for_ASE)
    .map{ left, right -> 
      def sampleId = left[0]
      def bam = left[1]
      def bai = left[2]
      def vcf = right[1]
      tuple(sampleId, vcf, bam, bai)  
    }
    .set { grouped_vcf_bam_bai_ch }
/* 
 * Process 6C: Allele-Specific Expression analysis with GATK ASEReadCounter.
 *             Calculates allele counts at a set of positions after applying 
 *             filters that are tuned for enabling allele-specific expression 
 *             (ASE) analysis
 */
process '6C_ASE_knownSNPs' {
  tag "$sampleId"
  publishDir "$params.results/$sampleId" 
  
  input:
      file genome from genome_file 
      file index from genome_index_ch
      file dict from genome_dict_ch
      set sampleId, file(vcf),  file(bam), file(bai) from grouped_vcf_bam_bai_ch
  
  output:
      file "ASE.tsv"
  
  script:
  """
  echo "${bam.join('\n')}" > bam.list
    
  $GATK -R ${genome} \
          -T ASEReadCounter \
          -o ASE.tsv \
          -I bam.list \
          -sites ${vcf}
  """
}
/*
 *  END OF PART 6
 ******/
