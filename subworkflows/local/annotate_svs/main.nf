include { ANNOTSV_ANNOTSV                 } from '../../../modules/nf-core/annotsv/annotsv/main'
include { SVDB_QUERY                      } from '../../../modules/nf-core/svdb/query/main'
include { ENSEMBLVEP_VEP as ENSEMBLVEP_SV } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow ANNOTATE_SVS {
    take:
    ch_vcf // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_fasta // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_sv_dbs // channel: [mandatory] [ val(meta), path(csv) ]
    ch_vep_cache // channel: [mandatory] [ val(meta), path(cache) ]
    val_vep_cache_version // string: [mandatory] default: 116
    ch_vep_extra_files // channel: [mandatory] [ path(files) ]
    val_annotate_annotsv // bool: [mandatory] should AnnotSV annotate family SVs
    ch_annotsv_annotations // channel: [optional] [ val(meta), path(annotations directory) ]

    main:
    ch_svdb_in = ch_sv_dbs.multiMap { filename, in_freq_info_key, in_allele_count_info_key, out_freq_info_key, out_allele_count_info_key ->
        vcf_dbs: filename
        in_frqs: in_freq_info_key
        in_occs: in_allele_count_info_key
        out_frqs: out_freq_info_key
        out_occs: out_allele_count_info_key
    }

    // Annotate with SVDB VCF "databases"
    SVDB_QUERY(
        ch_vcf,
        ch_svdb_in.in_occs.toList(),
        ch_svdb_in.in_frqs.toList(),
        ch_svdb_in.out_occs.toList(),
        ch_svdb_in.out_frqs.toList(),
        ch_svdb_in.vcf_dbs.toList(),
        [],
    )

    ENSEMBLVEP_SV(
        SVDB_QUERY.out.vcf.map { meta, vcf -> [meta, vcf, []] },
        "GRCh38",
        "homo_sapiens",
        val_vep_cache_version,
        ch_vep_cache,
        ch_fasta,
        ch_vep_extra_files,
        [[], []],
    )

    ch_annotsv_tsv = channel.empty()
    if (val_annotate_annotsv) {
        // Report on family SVs before the clinical and research sets are split.
        ANNOTSV_ANNOTSV(
            ENSEMBLVEP_SV.out.vcf.join(ENSEMBLVEP_SV.out.tbi, failOnMismatch: true, failOnDuplicate: true).map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
            ch_annotsv_annotations,
            [[], []],
            [[], []],
            [[], []],
        )
        ch_annotsv_tsv = ANNOTSV_ANNOTSV.out.tsv.mix(ANNOTSV_ANNOTSV.out.unannotated_tsv)
    }

    emit:
    vcf         = ENSEMBLVEP_SV.out.vcf // channel: [ val(meta), path(vcf) ]
    tbi         = ENSEMBLVEP_SV.out.tbi // channel: [ val(meta), path(tbi) ]
    annotsv_tsv = ch_annotsv_tsv // channel: [ val(meta), path(tsv) ]
}
