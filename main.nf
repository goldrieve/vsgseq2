nextflow.enable.dsl = 2

params.cores = "4"
params.trinitymem = "20"
params.outdir = "results"
params.help = ""
params.requestedcpus = 4
params.samplesheet = "samplesheet.csv"

include { TRIM } from './modules/trim'
include { ASSEMBLE } from './modules/assemble'

ch_samplesheet = Channel.fromPath(params.samplesheet)

ch_reads = ch_samplesheet.splitCsv(header:true).map {

    r1 = it['r1']
    r2 = it['r2']

    is_singleEnd = r2.toString()=='' ? true : false

    meta = [id: it['sample'], single_end: is_singleEnd]

    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

}

workflow {
   trimmed_reads = TRIM(ch_reads, params.cores)
   ASSEMBLE(ch_reads, params.cores, params.trinitymem)
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone!" : "Oops .. something went wrong")
}
