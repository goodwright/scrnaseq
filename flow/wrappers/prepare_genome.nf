// Test this workflow using:
// nextflow run ./flow/wrappers/prepare_genome.nf -profile docker,test -c ./nextflow.config -c ./flow/conf/test_wrapper.config -c ./flow/conf/prepare_genome.config --outdir ./results

include { PREPARE_SCRNASEQ } from '../../subworkflows/goodwright/prepare_genome/prepare_scrnaseq/main'

workflow {
    // Check files
    ch_fasta = file(params.fasta, checkIfExists: true)
    ch_gtf   = file(params.gtf, checkIfExists: true)

    // Run wrapper
    PREPARE_SCRNASEQ (
        ch_fasta,
        ch_gtf,
        "standard" // kb_workflow
    )
}
