# sv_somatic_cns
Consensus calling of Somatic Structural variants from Paired WGS

## Description
Pipeline using multiple SV callers for consensus structural variant calling from tumor/normal sequencing data. 

## Usage
  ```
  # Run the whole pipeline and consensus of calls
  nextflow run iarcbioinfo/sv_somatic_cns-nf -r v1.0 \
  -profile singularity  --tn_file tn_pairs..txt \
  --input_folder $PWD/CRAM \
  --ref hs38DH.fa \
  --all_sv_cns \
  --output_folder results

  #Run Delly and manta only
  nextflow run iarcbioinfo/sv_somatic_cns-nf -r v1.0 \
  -profile singularity  --tn_file tn_pairs..txt \
  --input_folder $PWD/CRAM \
  --ref hs38DH.fa \
  --delly \
  --manta \
  --output_folder results_delly_manta
  ```

## Dependencies

1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2. External software:
	- [DELLY v0.8.7](https://github.com/dellytools/delly)
	- [SVABA v1.1.0](https://github.com/walaj/svaba)
	- [MANTA v1.6.0](https://github.com/Illumina/manta)
	- [SURVIVOR v1.0.7](https://github.com/fritzsedlazeck/SURVIVOR)
	- [BCFTOOLS v1.10](https://github.com/samtools/bcftools)
	- [SAMTOOLS v1.10](https://github.com/samtools/samtools)

You can avoid installing all the external software by only installing Docker or singularity.
See the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository for more information.


## Input (mandatory)

  | Type      | Description   |
  |-----------|---------------|
  | --input_folder    | Folder containing all BAM/CRAM files |  
  | --tn_file    | File containing the list of names of BAM files to be processed |
  |--ref         |  Fasta file of reference genome [hg38.fa], should be indexed [hg38.fa.fai]|
  Flags to run each SV caller combinations
  | --delly  | run the Delly SV caller |
  | --manta  | run the Manta SV caller |
  | --svaba  | run the SVaba SV caller |
  Short-cut to enable all sv callers plus consensus with survivor
  | --all_sv_cns | run Delly, Manta, SVaba and integration with SURVIVOR|


### Example of Tumor/Normal pairs file (--tn_file)
A text file tabular separated, with the following header:

```
sampleID	tumor	normal
sample1_T1	sample1_T.cram	sample1_N.cram
sample2_T1	sample2_T.cram	sample2_N.cram
sample3_T1	sample3_T.cram	sample3_N.cram
```

### Optional parameters

| Name      | type | Description     |
|-----------|---------------|-----------------|
|      --bam     |       [flag] |active bam mode [def:cram]|
|     --output_folder |  [string] |name of output folder |
|      --cpu          |[Integer] | Number of CPUs[def:2] |
|      --mem |        [Integer] | Max memory [def:16Gb] |  

## Output
```
results
├── DELLY                               # DELLY result directory
│   ├── ...
├── SVABA                               # SVABA result directory
│   ├── ...
├── MANTA                               # MANTA result directory
│   ├── ...
├── SURVIVOR                            # SURVIVOR result directory
│   ├── ...
├── nf-pipeline_info                   # NEXTFLOW logs
```



## Limitations

The current version of the pipeline can handle Tumor/Normal pairs only but not multi-region WGS data.

## Common errors

### Singularity
The first time that the container is built from the docker image, the TMPDIR  should be defined in a non parallel file-system, you can set this like:

```
export TMPDIR=/tmp
```

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------|
  | Matthieu Foll*    |            follm@iarc.fr | Developer to contact for support (link to specific gitter chatroom) |
  | Alex Di Genova | digenovaa@fellows.iarc.fr| Developer |
