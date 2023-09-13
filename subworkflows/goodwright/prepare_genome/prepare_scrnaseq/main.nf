//
// Prepares genome reference files for nf-core/scrnaseq
//

include { GUNZIP as GUNZIP_FASTA } from '../../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_GTF   } from '../../../../modules/nf-core/gunzip/main.nf'
include { CELLRANGER_MKGTF       } from '../../../../modules/nf-core/cellranger/mkgtf/main'
include { CELLRANGER_MKREF       } from '../../../../modules/nf-core/cellranger/mkref/main'
include { KALLISTOBUSTOOLS_REF   } from '../../../../modules/nf-core/kallistobustools/ref/main'
include { STAR_GENOMEGENERATE    } from '../../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_SCRNASEQ {
    take:
    fasta       // file: /path/to/genome.fasta
    gtf         // file: /path/to/genome.gtf
    kb_workflow // string: should default to standard

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.toString().endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( [ [id:fasta.baseName], fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.of([ [id:fasta.baseName], fasta ])
    }
    //ch_fasta | view

    //
    // MODULE: Uncompress genome gtf file if required
    //
    ch_gtf = Channel.empty()
    if (gtf.toString().endsWith(".gz")) {
        ch_gtf      = GUNZIP_GTF ( [ [id:gtf.baseName], gtf ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.of([ [id:gtf.baseName], gtf ])
    }
    //ch_gtf | view

    //
    // MODULE: Create Cell ranger compatable GTF file - Filter GTF based on gene biotypes passed in params.modules
    //
    CELLRANGER_MKGTF(
        ch_gtf.collect{ it[1] }
    )
    ch_versions = ch_versions.mix(CELLRANGER_MKGTF.out.versions)

    //
    // MODULE: Make cell range reference genome
    //
    CELLRANGER_MKREF(
        ch_fasta.collect{ it[1] },
        CELLRANGER_MKGTF.out.gtf,
        "cellranger_reference"
    )
    ch_versions         = ch_versions.mix(CELLRANGER_MKREF.out.versions)
    ch_cellranger_index = CELLRANGER_MKREF.out.reference
    //ch_cellranger_index | view

    //
    // MODULE: Make bustools ref
    //
    KALLISTOBUSTOOLS_REF (
        ch_fasta.collect{ it[1] },
        ch_gtf.collect{ it[1] },
        kb_workflow
    )
    ch_versions       = ch_versions.mix(KALLISTOBUSTOOLS_REF.out.versions)
    ch_kallisto_index = KALLISTOBUSTOOLS_REF.out.index.collect()
    //ch_kallisto_index | view

    //
    // MODULE: Build star index
    //
    STAR_GENOMEGENERATE(
        ch_fasta.collect{ it[1] },
        ch_gtf.collect{ it[1] }
    )
    ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    ch_star_index = STAR_GENOMEGENERATE.out.index.collect()

    //
    // MODULE: Build salmon index
    //
    //TODO: SIMPLEAF is not on nf-core/modules yet

    emit:
    cellranger_index = ch_cellranger_index // channel: [ val(meta), [ index ] ]
    kallisto_index   = ch_kallisto_index   // channel: [ val(meta), [ index ] ]
    star_index       = ch_star_index       // channel: [ val(meta), [ index ] ]
    versions         = ch_versions         // channel: [ versions.yml ]
}
