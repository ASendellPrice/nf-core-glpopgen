//
// Population structure across all samples from joint genotype likelihoods
//

include { ANGSD_GL as ANGSD_GL_ALL } from '../../../modules/nf-core/angsd/gl/main'
include { PCANGSD_PCA              } from '../../../modules/local/pcangsd/pca/main'
include { ANGSD_NGSADMIX           } from '../../../modules/local/angsd/ngsadmix/main'

workflow BAM_TO_STRUCTURE {

    take:
    ch_bams       // channel: [ val(meta), [ bams ], [ bais ] ] single entry holding all samples
    ch_reference  // channel: [ fasta, fai ]
    admix_k_min   // integer: smallest K for NGSadmix
    admix_k_max   // integer: largest K for NGSadmix
    admix_reps    // integer: independent NGSadmix runs (seeds) per K

    main:

    if (admix_k_min < 1 || admix_k_max < admix_k_min) {
        error("BAM_TO_STRUCTURE: admix_k_min (${admix_k_min}) must be >= 1 and <= admix_k_max (${admix_k_max}).")
    }
    if (admix_reps < 1) {
        error("BAM_TO_STRUCTURE: admix_reps (${admix_reps}) must be >= 1.")
    }

    //
    // MODULE: Joint genotype likelihoods (Beagle) across all samples
    //
    ANGSD_GL_ALL (
        ch_bams,
        ch_reference.map { fasta, fai -> [ [:], fasta, fai ] },
        [ [:], [] ]     // error_file: only needed for -GL 4
    )

    //
    // MODULE: Genetic covariance matrix for PCA
    //
    PCANGSD_PCA (
        ANGSD_GL_ALL.out.genotype_likelihood
    )

    //
    // MODULE: Admixture proportions, one task per K and replicate (rep is used as the seed)
    //
    def ch_admix_input = ANGSD_GL_ALL.out.genotype_likelihood
        .flatMap { meta, beagle ->
            (admix_k_min..admix_k_max).collectMany { k ->
                (1..admix_reps).collect { rep -> [ meta + [ k: k, rep: rep ], beagle, k ] }
            }
        }

    ANGSD_NGSADMIX (
        ch_admix_input
    )

    emit:
    beagle     = ANGSD_GL_ALL.out.genotype_likelihood  // channel: [ val(meta), beagle.gz ]
    covariance = PCANGSD_PCA.out.covariance            // channel: [ val(meta), cov ]
    admix_qopt = ANGSD_NGSADMIX.out.qopt               // channel: [ val(meta + [ k, rep ]), qopt ]
    admix_log  = ANGSD_NGSADMIX.out.log                // channel: [ val(meta + [ k, rep ]), log ]
}
