// This ensures DSL2 syntax and process imports
nextflow.enable.dsl=2

include { DOWNLOAD_ERCC } from './modules/download.nf'
include { CONCAT_ERCC;
          SUBSAMPLE_GENOME;
          TO_PRED;
          TO_BED } from './modules/genome.nf'


include { DOWNLOAD_GUNZIP_REFERENCES } from './modules/download.nf'

def get_reference_urls(csv_table, organism_sci) {
    def organisms = [:]
      csv_table.splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    target_row = organisms[organism_sci.capitalize().replace("_"," ")]
    return tuple(organism_sci, target_row[2], target_row[3])
}

def get_annotation_url(csv_table, organism_sci) {
    def organisms = [:]
      csv_table.splitEachLine(",") {fields ->
          organisms[fields[1]] = fields
    }
    target_row = organisms[organism_sci.capitalize().replace("_"," ")]
    return target_row[6]
}

/**************************************************
* ACTUAL WORKFLOW  ********************************
**************************************************/
workflow references{
  take:
    organism_sci
    has_ercc
  main:
      if (params.ref_fasta && params.ref_gtf) {
        genome_annotations_pre_subsample = Channel.fromPath([params.ref_fasta, params.ref_gtf], checkIfExists: true).toList()
        genome_annotations_pre_subsample | view
      } else {
        // use assets table to find current fasta and gtf urls
        organism_sci | map{ org -> get_reference_urls( file(params.reference_table), org ) } | DOWNLOAD_GUNZIP_REFERENCES
        DOWNLOAD_GUNZIP_REFERENCES.out | set{ genome_annotations_pre_subsample }
      }
      // use assets table to find current annotations file
      organism_sci | map{ org -> get_annotation_url( file(params.reference_table), org ) } | set{ ch_gene_annotations_url }

      // SUBSAMPLING STEP : USED FOR DEBUG/TEST RUNS
      if ( params.genomeSubsample ) {
        SUBSAMPLE_GENOME( genome_annotations_pre_subsample, organism_sci )
        SUBSAMPLE_GENOME.out.build | set { genome_annotations_pre_ercc }
      } else {
        genome_annotations_pre_subsample | set { genome_annotations_pre_ercc }
      }

      // ERCC STEP : ADD ERCC Fasta and GTF to genome files
      DOWNLOAD_ERCC(has_ercc).ifEmpty([file("ERCC92.fa"), file("ERCC92.gtf")]) | set { ch_maybe_ercc_refs }
      CONCAT_ERCC( genome_annotations_pre_ercc, ch_maybe_ercc_refs, organism_sci, has_ercc )
      .ifEmpty { genome_annotations_pre_ercc.value }  | set { genome_annotations }


      TO_PRED( genome_annotations | map { it[1] }, organism_sci )
      TO_BED( TO_PRED.out, organism_sci )

  emit:
      genome_annotations = genome_annotations
      genome_bed = TO_BED.out
      gene_annotations = ch_gene_annotations_url
}
