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

        SVaba:
        --svaba_dbsnp                FILE        dbSNP file available at: https://data.broadinstitute.org/snowman/dbsnp_indel.vcf
        --svaba_targets              FILE        bed file with target positions for svaba
        --svaba_options              STRING      List of options to pass to svaba


    """.stripIndent()
}


// Show help message
if (params.help) exit 0, show_help()

//we load the tn_file for processing
if(params.tn_file == null) exit 0, show_help()
if(params.bam){
  params.ext=".bai"
}

//we create the channel for svaba, delly and manta
Channel.fromPath(returnFile(params.tn_file)).splitCsv(header: true, sep: '\t', strip: true)
                .map{row -> [ row.sampleID,
                              returnFile(params.input_folder + "/" +row.tumor),
                              returnFile(params.input_folder + "/" +row.tumor+params.ext),
                              returnFile(params.input_folder + "/" +row.normal),
                              returnFile(params.input_folder + "/" +row.normal+params.ext)]}
                  .into{genomes_svaba; genomes_delly; genomes_manta;}


running_tools = []
//we load the index file
fasta_ref = returnFile(params.ref)
fasta_ref_fai = returnFile( params.ref+'.fai' )

delly_blacklist = ""
manta_callable = ""

if (params.delly) {
    running_tools.add("Delly")
    delly_blacklist = returnFile("$baseDir/blacklist/human.hg38.excl.delly.tsv")
}

if (params.manta) {
    running_tools.add("Manta")
    manta_callable = returnFile("$baseDir/blacklist/manta_callable_chrs.hg38.bed.gz")
}


//bwa index for SVaba
fasta_ref_sa = null
fasta_ref_bwt = null
fasta_ref_ann = null
fasta_ref_amb = null
fasta_ref_pac = null
fasta_ref_alt = null
//bwa_index=Channel.create()
if (params.svaba) {
    running_tools.add("SVaba")
    //BWA index for SVaba
    //Channel.fromList(["BWAindex",returnFile( params.ref+'.sa'),
    //                  returnFile( params.ref+'.bwt' ),
    //                  returnFile( params.ref+'.ann' ),
    //                  returnFile( params.ref+'.amb' ),
    //                  returnFile( params.ref+'.pac' ),
    //                  returnFile( params.ref+'.alt' ),
    //                  ])
    //      .set{bwa_index;}

   fasta_ref_sa = returnFile( params.ref+'.sa' )
    fasta_ref_bwt = returnFile( params.ref+'.bwt' )
    fasta_ref_ann = returnFile( params.ref+'.ann' )
    fasta_ref_amb = returnFile( params.ref+'.amb' )
    fasta_ref_pac = returnFile( params.ref+'.pac' )
    fasta_ref_alt = returnFile( params.ref+'.alt' )
}


//aux funtions to check is a file exists
def returnFile(it) {
  // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "The following file: ${inputFile},  do not exist!!! see --help for more information"
    return inputFile
}

//delly process
process delly {
  cpus params.cpu
  memory params.mem+'G'
  tag {"Delly"+sampleID }

  publishDir "${params.output_folder}/DELLY/", mode: 'copy'

  input:
  set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from genomes_delly
  file fasta_ref
  file fasta_ref_fai
  file delly_blacklist

  output:
   //primary vcf file
   set val(sampleID), file("${sampleID}.delly_somatic.vcf") into delly_vcf
   //optional files
   file("*.{tsv,bcf}") into delly_output
  when: params.delly

  script:
  //run delly with mathched data
  """
  delly call -x  ${delly_blacklist}  -g ${fasta_ref} -o ${sampleID}.matched.bcf ${tumorBam} ${normalBam}
  #we call use the file
  bcftools view ${sampleID}.matched.bcf | grep "^#CHR" | awk '{print \$10"\ttumor"; print \$11"\tcontrol"}' > ${sampleID}.sample.tsv
  #we apply the somatic filter and keep only PASS variants
  delly filter -f somatic -s ${sampleID}.sample.tsv -p -o ${sampleID}.somatic.bcf ${sampleID}.matched.bcf
  #we convert the bcf file to VCF
  bcftools view ${sampleID}.somatic.bcf > ${sampleID}.delly_somatic.vcf
  """

}



//we run the SVaba caller
process svaba {
	cpus params.cpu
     memory params.mem+'G'
     tag { "SVABA"+sampleID }

     publishDir "${params.output_folder}/SVABA/", mode: 'copy'

     input :
     set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from genomes_svaba
     file fasta_ref
     file fasta_ref_fai
     //set val(sampleID),file(fasta_ref_sa),file(fasta_ref_bwt),file(fasta_ref_ann),file(fasta_ref_ann),file(fasta_ref_amb),file(fasta_ref_pac),file(fasta_ref_alt) from bwa_index
     file fasta_ref_sa
     file fasta_ref_bwt
     file fasta_ref_ann
     file fasta_ref_amb
     file fasta_ref_pac
     file fasta_ref_alt

     output:
     set val(sampleID), file("${sampleID}*.vcf") into svaba_vcf
     file "${sampleID}.alignments.txt.gz" into svaba_alignments

     when: params.svaba


     shell :
     if(params.targets) targets="-k ${params.svaba_targets}"
     else targets=""
     if(normalBam.baseName == 'None' ) normal=""
     else  normal="-n ${normalBam}"
     '''
     svaba run -t !{tumorBam} !{normal} -p !{params.cpu} !{params.svaba_dbsnp} -a somatic_run -G !{fasta_ref} !{targets} !{params.svaba_options}
     mv somatic_run.alignments.txt.gz !{sampleID}.alignments.txt.gz
     for f in `ls *.vcf`; do mv $f !{sampleID}.$f; done
     '''
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