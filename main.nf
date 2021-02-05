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
      --tumor_cram_folder   FOLDER                  Folder containing tumor CRAM files to be called.
      --normal_cram_folder  FOLDER                  Folder containing matched normal CRAM files.'

      --output_folder      FOLDER                   Output Folder [def:results]
    References
      --ref [file]                    Path to fasta reference
      --bwa_index [file]               Path to BWA index for running SVaba caller

      Profiles:

      -profile [str]              Configuration profile to use.
                                  Available: singularity

    Tool flags:
      --delly [bool]                  Run DELLY SV caller
      --manta [bool]                  Run Manta SV caller
      --SVaba [bool]                  Run SVaba SV caller

    Tool options:
        Delly:

        Manta:

        SVaba:
        --svaba_dbsnp                FILE        dbSNP file available at: https://data.broadinstitute.org/snowman/dbsnp_indel.vcf
        --svaba_targets              FILE        bed file with target positions for svaba
        --svaba_options              STRING      List of options to pass to svaba


    """.stripIndent()
}

//we init some params
params.help = null
params.cpu = 1
params.mem = 4


// Show help message
if (params.help) exit 0, show_help()

running_tools = []
if (params.delly) {
    running_tools.add("Delly")
    //reference.arriba = Channel.value(file(params.arriba_ref)).ifEmpty{exit 1, "Arriba reference directory not found!"}
}

if (params.manta) {
    running_tools.add("Manta")
    //reference.arriba = Channel.value(file(params.arriba_ref)).ifEmpty{exit 1, "Arriba reference directory not found!"}
}

if (params.svaba) {
    running_tools.add("SVaba")

    fasta_ref = returnFile(params.ref)
    //we ask if the BWA index exit
    fasta_ref_fai = returnFile( params.ref+'.fai' )
    fasta_ref_sa = returnFile( params.ref+'.sa' )
    fasta_ref_bwt = returnFile( params.ref+'.bwt' )
    fasta_ref_ann = returnFile( params.ref+'.ann' )
    fasta_ref_amb = returnFile( params.ref+'.amb' )
    fasta_ref_pac = returnFile( params.ref+'.pac' )
    fasta_ref_alt = returnFile( params.ref+'.alt' )

    //reference.arriba = Channel.value(file(params.arriba_ref)).ifEmpty{exit 1, "Arriba reference directory not found!"}
}






//aux funtions to check is a file exists
def returnFile(it) {
  // Return file if it exists
    inputFile = file(it)
    if (!file(inputFile).exists()) exit 1, "The following file: ${inputFile},  do not exist!!! see --help for more information"
    return inputFile
}


//we run the SVaba caller

process svaba {
	cpus params.cpu
     memory params.mem+'G'
     tag { sampleID }

     publishDir params.output_folder, mode: 'copy'

     input :
     set val(sampleID),file(tumorBam),file(tumorBai),file(normalBam),file(normalBai) from bams
     file fasta_ref
     file fasta_ref_fai
     file fasta_ref_sa
     file fasta_ref_bwt
     file fasta_ref_ann
     file fasta_ref_amb
     file fasta_ref_pac
     file fasta_ref_alt

     output:
     set val(sampleID), file("${sampleID}*.vcf") into vcf
     file "${sampleID}.alignments.txt.gz" into alignments

     shell :
     if(params.targets) targets="-k ${params.targets}"
     else targets=""
     if(normalBam.baseName == 'None' ) normal=""
     else  normal="-n ${normalBam}"
     '''
     svaba run -t !{tumorBam} !{normal} -p !{params.cpu} !{dbsnp_par} !{params.dbsnp} -a somatic_run -G !{fasta_ref} !{targets} !{params.options}
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
        F\u001b[31;1mU\u001b[32;1mS\u001b[33;1mI\u001b[0mO\u001b[33;1mN\u001b[31;1m : Gene\u001b[32;1m Fusion\u001b[33;1m Caller\u001b[31;1m (${workflow.manifest.version})
        """
}
