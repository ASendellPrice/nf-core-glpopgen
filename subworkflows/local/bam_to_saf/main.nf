//
// Estimate site allele frequency likelihoods (SAF) per population with ANGSD
//

include { ANGSD_SOAPSNPCALIBRATION } from '../../../modules/nf-core/angsd/soapsnpcalibration/main'
include { BEAGLE_SUBSET            } from '../../../modules/local/beagle/subset/main'
include { PCANGSD_INBREEDING       } from '../../../modules/nf-core/pcangsd/inbreeding/main'
include { ANGSD_DOSAF              } from '../../../modules/nf-core/angsd/dosaf/main'

workflow BAM_TO_SAF {

    take:
    ch_bams       // channel: [ val(meta), [ bams ], [ bais ] ] one entry per population
    ch_beagle     // channel: [ val(meta), beagle.gz ] all-samples genotype likelihoods; meta.samples gives the column order
    ch_reference  // channel: [ fasta, fai ]
    ch_ancestral  // channel: [ ancestral_fasta, ancestral_fai ] or [ [], [] ]
    gl_model      // integer: ANGSD -GL model (1, 2 or 3)
    saf_method    // integer: ANGSD -doSAF mode (1 = standard, 2 = inbreeding-aware)

    main:

    if (!(gl_model in [1, 2, 3])) {
        error("BAM_TO_SAF: unsupported gl_model '${gl_model}'. Use 1, 2 or 3 (-GL 4 requires an error file, which is not yet supported).")
    }
    if (!(saf_method in [1, 2])) {
        error("BAM_TO_SAF: unsupported saf_method '${saf_method}'. Use 1 or 2.")
    }

    //
    // MODULE: SOAPsnp calibration matrix per population (-GL 3 only)
    //
    def ch_calibration = ch_bams.map { meta, _bams, _bais -> [ meta, [] ] }
    if (gl_model == 3) {
        ANGSD_SOAPSNPCALIBRATION (
            ch_bams,
            ch_reference
        )
        ch_calibration = ANGSD_SOAPSNPCALIBRATION.out.soapsnp_calibration
    }

    //
    // MODULES: Per-sample inbreeding coefficients per population (-doSAF 2 only)
    // Each population's columns are taken from the all-samples Beagle file, in the same
    // sample order as its BAM list, so the -indF values line up with doSAF's bamlist
    //
    def ch_inbreeding     = ch_bams.map { meta, _bams, _bais -> [ meta, [] ] }
    def ch_beagle_samples = channel.empty()
    if (saf_method == 2) {
        def ch_subset_input = ch_bams
            .combine(ch_beagle)
            .map { meta, _bams, _bais, meta_all, beagle -> [ meta, beagle, meta_all.samples, meta.samples ] }

        BEAGLE_SUBSET (
            ch_subset_input
        )
        PCANGSD_INBREEDING (
            BEAGLE_SUBSET.out.beagle
        )
        ch_inbreeding     = PCANGSD_INBREEDING.out.inbreeding_coefficients
        ch_beagle_samples = BEAGLE_SUBSET.out.samples
    }

    //
    // Join optional per-population inputs by meta so each stays paired with its own BAMs
    //
    def ch_dosaf_input = ch_bams
        .join(ch_calibration, failOnMismatch: true)
        .join(ch_inbreeding, failOnMismatch: true)
        .multiMap { meta, bams, bais, calibration, inbreeding ->
            bams:        [ meta, bams, bais ]
            calibration: calibration
            inbreeding:  inbreeding
        }

    //
    // MODULE: Site allele frequency likelihoods per population
    //
    ANGSD_DOSAF (
        ch_dosaf_input.bams,
        ch_reference,
        ch_ancestral,
        [],                 // error_file: only needed for -GL 4
        ch_dosaf_input.inbreeding,
        ch_dosaf_input.calibration,
        gl_model,
        saf_method
    )

    emit:
    saf                     = ANGSD_DOSAF.out.saf                   // channel: [ val(meta), saf.idx, saf.pos.gz, saf.gz ]
    inbreeding_coefficients = ch_inbreeding                         // channel: [ val(meta), inbreed.samples ] or [ val(meta), [] ]
    beagle_samples          = ch_beagle_samples                     // channel: [ val(meta), samples.txt ] (-doSAF 2 only)
}
