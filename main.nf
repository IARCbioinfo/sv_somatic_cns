#!/usr/bin/env nextflow


//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarc/sv_somatic_cns --tn_file pairs.txt -profile singularity

    Mandatory arguments:
      --tn_file  [file]                             file with tabular data for each sample to process [sampleID tumor normal]
      --input_folder   FOLDER                       Folder containing CRAM/BAM files to be called.
      --output_folder      FOLDER                   Output Folder [def:results]
      --bam                                         File to process are BAM [def:CRAM]
    References
      --ref [file]                    Path to fasta reference
      --bwa_index [file]               Path to BWA index for running SVaba caller

      Profiles:

      -profile [str]              Configuration profile to use.
                                  Available: singularity

    Tool flags:
      --delly [bool]                  Run DELLY SV caller
      --manta [bool]                  Run Manta SV caller
      --svaba [bool]                  Run SVaba SV caller

    Tool options:
        Delly:

        Manta:
        --manta_cfg                   Manta config file [def:"/opt/conda/envs/sv_somatic_cns/share/manta-1.6.0-0/bin/configManta.py"]
        --manta_samtools              path to samtools [def:"/opt/conda/envs/sv_somatic_cns/bin/samtools"]
        SVaba:
        --svaba_dbsnp                FILE        dbSNP file available at: https://data.broadinstitute.org/snowman/dbsnp_indel.vcf
        --svaba_targets              FILE        bed file with target positions for svaba
        --svaba_options              STRING      List of options to pass to svaba
        --svaba_by_chr  [bool]        Run SVABA by Chromsosome [def:true]           

    """.stripIndent()
}


// Show help message
if (params.help) exit 0, show_help()

//we load the tn_file for processing
if(params.tn_file == null) exit 0, show_help()
if(params.ref == null) exit 0, show_help()
if(params.bam){
  params.ext=".bai"
}

//we create the channel for svaba, delly and manta
Channel.fromPath(returnFile(params.tn_file)).splitCsv(header: true, sep: '\t', strip: true)
                .map{row -> [ row.sampleID,
                              file(params.input_folder + "/" +row.tumor),
                              file(params.input_folder + "/" +row.tumor+params.ext),
                              file(params.input_folder + "/" +row.normal),
                              file(params.input_folder + "/" +row.normal+params.ext)]}
                .into{genomes_svaba; genomes_delly; genomes_manta}


running_tools = []
fasta_ref_c = Channel.value(returnFile(params.ref)).ifEmpty{exit 1, "Fasta file not found: ${params.ref}"}
fasta_ref_fai_c = Channel.value(returnFile( params.ref+'.fai' )).ifEmpty{exit 1, "FAI file not found: ${params.ref}.fai"}



//channels for callabe or blacklist regions hg38
manta_callable_c = false
delly_blacklist_c = false
manta_callable_i = false
svaba_blacklist_c = false

if (params.delly) {
    running_tools.add("Delly")
    delly_blacklist_c = Channel.value(returnFile("$baseDir/blacklist/human.hg38.excl.delly.tsv")).ifEmpty{exit 1, "File not found: $baseDir/blacklist/human.hg38.excl.delly.tsv"}
}

if (params.manta) {
    running_tools.add("Manta")
    manta_callable_c = Channel.value(returnFile("$baseDir/blacklist/manta_callable_chrs.hg38.bed.gz")).ifEmpty{exit 1, "File not found: $baseDir/blacklist/manta_callable_chrs.hg38.bed.gz"}
    manta_callable_i = Channel.value(returnFile("$baseDir/blacklist/manta_callable_chrs.hg38.bed.gz.tbi")).ifEmpty{exit 1, "File not found: $baseDir/blacklist/manta_callable_chrs.hg38.bed.gz.tbi"}
}


//bwa index for SVaba
fasta_ref_sa = ""
fasta_ref_bwt = ""
fasta_ref_ann = ""
fasta_ref_amb = ""
fasta_ref_pac = ""
fasta_ref_alt = ""
//bwa_index=Channel.create()
if (params.svaba) {
    running_tools.add("SVaba")
    //BWA index for SVaba
   fasta_ref_sa = returnFile( params.ref+'.sa' )
   fasta_ref_bwt = returnFile( params.ref+'.bwt' )
   fasta_ref_ann = returnFile( params.ref+'.ann' )
   fasta_ref_amb = returnFile( params.ref+'.amb' )
   fasta_ref_pac = returnFile( params.ref+'.pac' )
   fasta_ref_alt = returnFile( params.ref+'.alt' )
   svaba_blacklist_c = Channel.value(returnFile("$baseDir/blacklist/human.hg38.excl.svaba.bed")).ifEmpty{exit 1, "File not found: $baseDir/blacklist/human.hg38.excl.svaba.bed"}
}



//manta process

process manta  {
  cpus params.cpu
  memory params.mem+'G'
  tag "${sampleID}-Manta"

  publishDir "${params.output_folder}/MANTA/", mode: 'copy'

  input:
  set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from genomes_manta
  file(fasta_ref) from fasta_ref_c
  file(fasta_ref_fai) from fasta_ref_fai_c
  file (manta_callable) from manta_callable_c
  file (manta_callable_index) from manta_callable_i

  output:
   //primary vcf file
   set val(sampleID), file("${sampleID}.manta_somatic_inv.vcf") into manta_vcf
   //set val(sampleID), file("${sampleID}.manta_somatic.vcf") into manta_vcf
   //additional files
   file("${sampleID}.manta_somatic.vcf") into manta_output
   file("${sampleID}_results") into manta_res_dir
  when: params.manta

  script:
  //we set the computational resources for this tool
  if (params.debug==false){
  """
  #we create the configuration file
  CFG_MANTA=${params.manta_cfg}
  python \$CFG_MANTA \\
  --normalBam=${normalBam} \\
  --tumorBam=${tumorBam} \\
  --runDir=${sampleID}.matched \\
  --referenceFasta=${fasta_ref} \\
  --callRegions=${manta_callable}
  #we run the manta caller
  python ${sampleID}.matched/runWorkflow.py  --quiet -j ${params.cpu} -g ${params.mem}
  #we recover the result files
  cp ${sampleID}.matched/results/variants/somaticSV.vcf.gz   ${sampleID}.manta_somatic.vcf.gz
  gzip -dc ${sampleID}.manta_somatic.vcf.gz > ${sampleID}.manta_somatic.vcf
  python ${baseDir}/aux_scripts/manta_convertINV.py ${params.manta_samtools} ${fasta_ref} ${sampleID}.manta_somatic.vcf > ${sampleID}.manta_somatic_inv.vcf
  mv ${sampleID}.matched/results ${sampleID}_results
  #mv ${sampleID}.matched/results/variants/candidateSmallIndels.vcf.gz ${sampleID}.manta_candidateSmallIndels.vcf.gz
  #mv ${sampleID}.matched/results/variants/candidateSV.vcf.gz ${sampleID}.manta_candidateSV.vcf.gz
  #mv ${sampleID}.matched/results/variants/diploidSV.vcf.gz ${sampleID}.manta_diploidSV.vcf.gz
  """
  }else{
    """
    touch ${sampleID}.manta_somatic_inv.vcf
    touch ${sampleID}.manta_somatic.vcf
    mkdir ${sampleID}_results
    """
  }
}


//delly process
process delly {
  cpus 1
  memory params.mem+'G'
  tag "${sampleID}-delly"

  publishDir "${params.output_folder}/DELLY/", mode: 'copy'

  input:
  set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from genomes_delly
  file(fasta_ref) from fasta_ref_c
  file(fasta_ref_fai) from fasta_ref_fai_c
  file (delly_blacklist) from delly_blacklist_c

  output:
   //primary vcf file
   set val(sampleID), file("${sampleID}.delly_somatic.vcf") into delly_vcf
   //optional files
   file("*.{tsv,bcf}") into delly_output
  when: params.delly

  script:
  //run delly with mathched data
  if (params.debug==false){
  """
  delly call -x  ${delly_blacklist}  -g ${fasta_ref} -o ${sampleID}.matched.bcf ${tumorBam} ${normalBam}
  #we call use the file
  bcftools view ${sampleID}.matched.bcf | grep "^#CHR" | awk '{print \$10"\ttumor"; print \$11"\tcontrol"}' > ${sampleID}.sample.tsv
  #we apply the somatic filter and keep only PASS variants
  delly filter -f somatic -s ${sampleID}.sample.tsv -p -o ${sampleID}.somatic.bcf ${sampleID}.matched.bcf
  #we convert the bcf file to VCF
  bcftools view ${sampleID}.somatic.bcf > ${sampleID}.delly_somatic.vcf
  """
  }else{
    """
      touch ${sampleID}.sample.tsv
      touch ${sampleID}.matched.bcf
      touch ${sampleID}.somatic.bcf
      touch ${sampleID}.delly_somatic.vcf
    """
  }
}

//we run the SVaba caller
process svaba {
	   cpus params.cpu
     memory params.mem+'G'
     tag "${sampleID}-SVaba"

     publishDir "${params.output_folder}/SVABA/", mode: 'copy'

     input:
     set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from genomes_svaba
     file(fasta_ref) from fasta_ref_c
     file(fasta_ref_fai) from fasta_ref_fai_c
     file fasta_ref_sa
     file fasta_ref_bwt
     file fasta_ref_ann
     file fasta_ref_amb
     file fasta_ref_pac
     file fasta_ref_alt
     file (svaba_blacklist) from svaba_blacklist_c

     output:
     set val(sampleID), file("${sampleID}*.vcf") into svaba_vcf
     file "${sampleID}.{alignments.txt.gz,bps.txt.gz,log}" into svaba_alignments

     when: params.svaba

     script:
     if(params.svaba_targets) targets="-k ${params.svaba_targets}"
     else targets=""
     if(normalBam.baseName == 'None' ) normal=""
     else  normal="-n ${normalBam}"
     if(params.svaba_dbsnp == "None") dbsnp=""
     else dbsnp="--dbsnp-vcf ${params.svaba_dbsnp}"
    if (params.debug==false){
      if(params.svaba_by_chr == true){
        //we run svaba by chromosome
        """
        #we define the chromsomes to call
        CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        for chr in \$CHRS ; do
        svaba run -k \$chr -t ${tumorBam} ${normal} -p ${params.cpu} ${dbsnp} -B ${svaba_blacklist} \\
        -a somatic_run_\$chr -G ${fasta_ref} ${targets} ${params.svaba_options} ;
        done
        #we merge the results into a single genome
        perl ${baseDir}/aux_scripts/merge_chrs_svaba.pl -s ${sampleID}
        #generates the output files
        """

      }else{
     """
     svaba run -t ${tumorBam} ${normal} -p ${params.cpu} ${dbsnp} -B ${svaba_blacklist} -a somatic_run -G ${fasta_ref} ${targets} ${params.svaba_options}
     mv somatic_run.alignments.txt.gz ${sampleID}.alignments.txt.gz
     mv somatic_run.bps.txt.gz ${sampleID}.bps.txt.gz
     mv somatic_run.log ${sampleID}.log
     for f in `ls *.vcf`; do mv \$f ${sampleID}.\$f; done
     """
     }
   }else{
     """
     touch ${sampleID}.alignments.txt.gz
     touch ${sampleID}.somatic.sv.vcf
     touch ${sampleID}.somatic.indel.vcf
     touch ${sampleID}.germline.sv.vcf
     touch ${sampleID}.germline.indel.vcf
     """
   }
}

//svaba run -t T680_DA_C000KIJ_H2M5HCCX2_hs38dh_MERGE_PE_3-4.reliable.cram -n T680_DA_C000KII_H2MKNCCX2_hs38dh_MERGE_PE_8.reliable.cram -p 2 --dbsnp-vcf /data/scratch/digenovaa/lungNENomics/SVs/databases/dbsnp_146.hg38.indels.vcf -B human.hg38.excl.svaba.bed -a s_run_chr1 -G hs38DH.fa -k chr2
//s_run_chr1.alignments.txt.gz  s_run_chr1.svaba.germline.indel.vcf	      s_run_chr1.svaba.unfiltered.germline.sv.vcf
//s_run_chr1.bps.txt.gz	      s_run_chr1.svaba.germline.sv.vcf		      s_run_chr1.svaba.unfiltered.somatic.indel.vcf
//s_run_chr1.contigs.bam	      s_run_chr1.svaba.somatic.indel.vcf	      s_run_chr1.svaba.unfiltered.somatic.sv.vcf
//s_run_chr1.discordant.txt.gz  s_run_chr1.svaba.somatic.sv.vcf
//s_run_chr1.log		      s_run_chr1.svaba.unfiltered.germline.indel.vcf

//aux funtions to check is a file exists
def returnFile(it) {
  // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "The following file: ${inputFile},  do not exist!!! see --help for more information"
      return inputFile
}





//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        Consensus Somatic Caller : (v${workflow.manifest.version})
        """
}
