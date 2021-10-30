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
      --ref [file]                    Path to fasta reference including BWA index

      Profiles:

      -profile [str]              Configuration profile to use.
                                  Available: singularity

    Tool flags:
      --delly [bool]                  Run DELLY SV caller
      --manta [bool]                  Run Manta SV caller
      --svaba [bool]                  Run SVaba SV caller
      --all_sv_cns [bool]             Run all SV callers plus consensus
    Tool options:
        Delly:

        Manta:
        --manta_cfg                   Manta config file [def:"/opt/conda/envs/sv_somatic_cns/share/manta-1.6.0-0/bin/configManta.py"]
        --manta_samtools              path to samtools [def:"/opt/conda/envs/sv_somatic_cns/bin/samtools"]
        SVaba:
        --svaba_dbsnp                FILE        dbSNP file available at: https://data.broadinstitute.org/snowman/dbsnp_indel.vcf
        --svaba_targets              FILE        bed file with target positions for svaba
        --svaba_options              STRING      List of options to pass to svaba
        --svaba_by_chr  [bool]        Run SVABA by Chromosome recommended when alignments are in CRAM format [def:true]
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

if (params.delly || params.all_sv_cns) {
    running_tools.add("Delly")
    delly_blacklist_c = Channel.value(returnFile("$baseDir/blacklist/human.hg38.excl.delly.tsv")).ifEmpty{exit 1, "File not found: $baseDir/blacklist/human.hg38.excl.delly.tsv"}
}

if (params.manta || params.all_sv_cns) {
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

if (params.svaba || params.all_sv_cns) {
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
   set val(sampleID), file("${sampleID}.manta_somatic_inv.pass.vcf") into manta_vcf
   //additional files
   file("${sampleID}.manta_somatic_inv.pass.vcf") into manta_output
   file("${sampleID}_results") into manta_res_dir
  when: params.manta || params.all_sv_cns

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
  bcftools view -f PASS -o ${sampleID}.manta_somatic_inv.pass.vcf -O v ${sampleID}.manta_somatic_inv.vcf
  mv ${sampleID}.matched/results ${sampleID}_results
  """
  }else{
    """
    touch ${sampleID}.manta_somatic_inv.vcf
    touch ${sampleID}.manta_somatic.vcf
    touch ${sampleID}.manta_somatic_inv.pass.vcf
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
  when: params.delly || params.all_sv_cns

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
     set val(sampleID), file("${sampleID}.svaba.somatic.sv.types.vcf") into svaba_vcf
     file "${sampleID}.{alignments.txt.gz,bps.txt.gz,log}" into svaba_alignments optional true
     file "${sampleID}.svaba.unfiltered*.vcf" into svaba_unfiltered
     file "${sampleID}.svaba.germline*.vcf" into svaba_germline
     file "${sampleID}.{svaba.somatic.indel.vcf,svaba.somatic.sv.vcf}" into svaba_somatic

     when: params.svaba || params.all_sv_cns

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
        #we define the chromsomes to call, here we assume hg38 reference
        CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
        for chr in \$CHRS ; do
        svaba run -k \$chr -t ${tumorBam} ${normal} -p ${params.cpu} ${dbsnp} -B ${svaba_blacklist} \\
        -a somatic_run_\$chr -G ${fasta_ref} ${targets} ${params.svaba_options} ;
        done
        #we merge the results into a single genome
        perl ${baseDir}/aux_scripts/merge_chrs_svaba.pl -s ${sampleID} -p somatic_run
        #generates the output files
        #we remove the duplicated translocations that results for runing svABA by chromosome
        perl ${baseDir}/aux_scripts/remove_duplicated_translocations_svaba_by_chr.pl -a  ${sampleID}.svaba.somatic.sv.vcf >  ${sampleID}.svaba.somatic.sv.dedup.vcf
        #we add the types to svABA predictions
        perl  ${baseDir}/aux_scripts/add_type_svaba.pl -a ${sampleID}.svaba.somatic.sv.dedup.vcf > ${sampleID}.svaba.somatic.sv.types.vcf
        """

      }else{
     """
     svaba run -t ${tumorBam} ${normal} -p ${params.cpu} ${dbsnp} -B ${svaba_blacklist} -a somatic_run -G ${fasta_ref} ${targets} ${params.svaba_options}
     mv somatic_run.alignments.txt.gz ${sampleID}.alignments.txt.gz
     mv somatic_run.bps.txt.gz ${sampleID}.bps.txt.gz
     mv somatic_run.log ${sampleID}.log
     mv
     for f in `ls *.vcf`; do mv \$f ${sampleID}.\${f/somatic_run.}; done

     perl  ${baseDir}/aux_scripts/add_type_svaba.pl -a ${sampleID}.svaba.somatic.sv.dedup.vcf > ${sampleID}.svaba.somatic.sv.types.vcf
     """
     }
   }else{
     """
     touch ${sampleID}.alignments.txt.gz
     touch ${sampleID}.svaba.somatic.sv.vcf
     touch ${sampleID}.svaba.somatic.indel.vcf
     touch ${sampleID}.svaba.germline.sv.vcf
     touch ${sampleID}.svaba.germline.indel.vcf
     touch ${sampleID}.svaba.somatic.sv.types.vcf
     touch ${sampleID}.svaba.unfiltered.germline.indel.vcf
     touch ${sampleID}.svaba.unfiltered.germline.sv.vcf
     touch ${sampleID}.svaba.unfiltered.somatic.indel.vcf
     touch ${sampleID}.svaba.unfiltered.somatic.sv.vcf
     """
   }
}


//we merge integrate the calls using VURVIVOR
//merge manta delly
m_m_d=delly_vcf.join(manta_vcf)
survivor_input=m_m_d.join(svaba_vcf)

//m_m_d_s.view()

process SURVIVOR{
  cpus params.cpu
  memory params.mem+'G'
  tag "${sampleID}-SURVIVOR"

  publishDir "${params.output_folder}/SURVIVOR/", mode: 'copy'
  input:
  set val(sampleID),file(delly_v),file(manta_v),file(svaba_v) from survivor_input
  output:
   file("*.cns.*") into survivor_output
  script:
  if (params.debug==false){
   """
   #we run SURVIVOR plus some filters
   perl ${baseDir}/aux_scripts/merge_callers_survivor_matched.pl -a ${manta_v} \\
   -b ${delly_v} -c ${svaba_v} -p ${sampleID}.cns
   #Veen data for consensus with at least two tools and at least 15 pair-end read support for single-tool predictions
   sh ${baseDir}/aux_scripts/get_data_veen.sh ${sampleID}.cns.integration.vcf ${sampleID}.cns.veen.integration.txt
   #Veen diagram for
   sh ${baseDir}/aux_scripts/get_data_veen.sh ${sampleID}.cns.survivor.vcf ${sampleID}.cns.veen.survivor.txt
   # We convert the *.cns.integration.vcf files to bedpe format
   SURVIVOR vcftobed ${sampleID}.cns.integration.vcf -1 -1 ${sampleID}.cns.integration.bedpe
   """
  }else{
    """
    #we run SURVIVOR plus some filters
    echo perl ${baseDir}/aux_scripts/merge_callers_survivor_matched.pl -a ${manta_v} \\
                -b ${delly_v} -c ${svaba_v} -p ${sampleID}.cns
    #Veen data for consensus with at least two tools and at least 15 pair-end read support for single-tool predictions
    echo sh ${baseDir}/aux_scripts/get_data_veen.sh ${sampleID}.cns.integration.vcf ${sampleID}.cns.veen.integration.txt
    #Veen diagram for
    echo sh ${baseDir}/aux_scripts/get_data_veen.sh ${sampleID}.cns.survivor.vcf ${sampleID}.cns.veen.survivor.txt
    # We convert the *.cns.integration.vcf files to bedpe format
    echo SURVIVOR vcftobed ${sampleID}.cns.integration.vcf -1 -1 ${sampleID}.cns.integration.bedpe
    touch ${sampleID}.cns.veen.integration.txt
    touch ${sampleID}.cns.veen.survivor.txt
    touch ${sampleID}.cns.integration.vcf
    touch ${sampleID}.cns.survivor.vcf
    touch ${sampleID}.cns.lst
    touch ${sampleID}.cns.survivor.log
    touch ${sampleID}.cns.integration.bedpe
    """
  }

}


//aux funtions to check is a file exists
def returnFile(it) {
  // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "The following file: ${inputFile},  do not exist!!! see --help for more information"
      return inputFile
}


// print the calling parameter to the log and a log file
def print_params () {
  //software versions for v2.0
 def software_versions = ['delly' : '0.8.7',
                          'manta' : '1.6.0',
                          'svaba' : '1.1.0',
                          'survivor' : '1.0.7',
                          'bcftools' : '1.10',
                          'samtools' : '1.10']
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
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
