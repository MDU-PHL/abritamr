#!/usr/bin/env nextflow
nextflow.preview.dsl=2
version = '1.0'


// set default value for the contigs channel
contigs = Channel.fromPath("*/contigs.fa")
                    .map { files -> tuple([id:files.parent.name, species_exp: "${params.species}"], files)}

// to run for mdu qc... You need to have a very good reason to do this!!
// if( params.mduqc == 'mduqc' ) {
//     def job_types = [:]
//     def species_exp = [:]
//     distribute_table = file("distribute_table.txt")
//     distribution_reader =   distribute_table.newReader()
//     // read the file
//     distribution_reader.eachLine { line ->
//         if( line.split('\t')[0] != 'MDUID') {
//             job_types[line.split('\t')[0]] = line.split('\t')[3]
//             species_exp[line.split('\t')[0]] = line.split('\t')[1]
//         }
//     }
//     contigs = Channel.fromPath("QC/*/contigs.fa")
//                     .map { files -> tuple([id:files.parent.name,job_type: 'standard',species_exp: species_exp[files.parent.name], files)}
// } else {
//     contigs = Channel.fromPath("*/contigs.fa")
//                     .map { files -> tuple([id:files.parent.name, species_exp: $params.species], files)}
// }

params.species_options = ["Acinetobacter_baumannii", "Campylobacter", "Enterococcus_faecalis", "Enterococcus_faecium", "Escherichia", "Klebsiella", "Salmonella", "Staphylococcus_aureus", "Staphylococcus_pseudintermedius", "Streptococcus_agalactiae", "Streptococcus_pneumoniae", "Streptococcus_pyogenes", "Vibrio_cholerae"]

workflow {

    include { AMRFINDER;COLLATE_AMRFINDER_ISOLATE;COLLATE_AMRFINDER_MATCHES;COLLATE_AMRFINDER_PARTIALS } from './modules/main' addParams( options: [args:"" ,args2: 1] )

    AMRFINDER( contigs )
    COLLATE_AMRFINDER_ISOLATE( AMRFINDER.out.amrfinder )

    abritamr_matches = COLLATE_AMRFINDER_ISOLATE.out.abritamr_matches_isolates
                                .map { cfg, amrfinder_matches-> amrfinder_matches }
                                .collect()
    abritamr_partials = COLLATE_AMRFINDER_ISOLATE.out.abritamr_partials_isolates                              
                                .map { cfg, amrfinder_partials -> amrfinder_partials }
                                .collect()
    COLLATE_AMRFINDER_MATCHES ( abritamr_matches )
    COLLATE_AMRFINDER_PARTIALS ( abritamr_partials )


}