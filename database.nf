#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.threads = 20
params.out = "${launchDir}/data"

workflow {
    // Define data sources
    def foodb = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
    def genbank_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    def taxdump = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    //
    // 1) Core DB downloads & matching
    //
    download(foodb, genbank_summary) | 
        get_taxids
    download_taxa_dbs(taxdump)
    get_lineage(get_taxids.out.combine(download_taxa_dbs.out)) | 
        match_taxids

    //
    // 2) Contig (nucleotide) download: one monolithic job
    //
    match_taxids.out | 
        download_sequences

    // Flatten the tuple of paths into a simple seq channel
    contig_seqs = download_sequences.out.map{ it[0] }.flatten()

    //
    // 3) Read matches.csv that was CREATED by match_taxids into a channel for genomes
    //
    matches_ch = match_taxids.out |
        splitCsv(header:true) |
        map { row ->
            [
                id:    row.id,
                db:    row.db,
                url:   row.url,
                taxid: row.matched_taxid.toInteger()
            ]
        }

    //
    // 4) Per‑assembly GenBank downloads
    //
    genbank_ch = matches_ch
        .filter { it.db == 'genbank' }
        .map    { hit -> tuple(hit.id, hit.url, hit.taxid) }

    download_genome(genbank_ch)

    // Collect genome fna.gz files
    genome_seqs = download_genome.out

    //
    // 5) Merge contigs + genomes and proceed
    //
    // `merge` will interleave; you can also use `concat` if ordering matters
    all_seqs = contig_seqs.merge(genome_seqs)

    // Merge manifest files  
    // Collect only genome manifest files for make_manifest 
    // Update the make_manifest call to use both manifest file collections
    make_manifest(
        download_sequences.out.map{ it[1] },
        download_genome.out.map{ it[1] }.collect()
    )
    
    all_seqs | sketch

    ani(sketch.out.collect())

    food_mappings(match_taxids.out)
}


process download {
    cpus 1
    publishDir "${params.out}/dbs", mode: 'copy'

    input:
    val foodb
    val genbank_summary

    output:
    tuple path("foodb"), path("genbank_summary.tsv")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${foodb} -O foodb.tgz && \\
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${genbank_summary} -O genbank_summary.tsv && \\
    tar -xf foodb.tgz && \\
    mv foodb_*_csv foodb
    """
}

process get_taxids {
    cpus 1

    input:
    tuple path(foodb), path(gb_summary)

    output:
    tuple path("foodb"), path("taxids.tsv"), path("${gb_summary}")

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    dt <- fread("${gb_summary}", sep="\t")[
        grepl("ftp.ncbi.nlm.nih.gov", ftp_path, fixed = TRUE)
    ]
    genbank <- dt[!is.na(taxid), .(taxid = as.character(unique(taxid)))]
    genbank[, "source" := "genbank"]
    dt <- fread("${foodb}/Food.csv")
    foodb <- dt[!is.na(ncbi_taxonomy_id), .(taxid = ncbi_taxonomy_id)]
    foodb[, "source" := "foodb"]
    fwrite(rbind(genbank, foodb), "taxids.tsv", col.names=F, sep="\t")
    """
}

process download_taxa_dbs {
    cpus 1

    input:
    val(taxdump)

    output:
    path("taxdump")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \\
        ${taxdump} && \\
        mkdir taxdump && tar -xf taxdump.tar.gz --directory taxdump
    """
}

process get_lineage {
    cpus 1

    input:
    tuple path(foodb), path(taxids), path(gb_summary), path(taxadb)

    output:
    tuple path("$foodb"), path("lineage.txt"), path("lineage_ids.txt"), path("${gb_summary}")

    script:
    """
    taxonkit lineage --data-dir $taxadb -i 1 $taxids > raw.txt && \\
    taxonkit reformat --data-dir $taxadb -i 3 raw.txt > lineage.txt && \\
    taxonkit reformat --data-dir $taxadb -t -i 3 raw.txt > lineage_ids.txt
    """
}

process match_taxids {
    cpus 1
    publishDir params.out, mode: 'copy'

    input:
    tuple path(foodb), path(lineage), path(lineage_ids), path(gb_summary)

    output:
    path("matches.csv")

    script:
    """
    match.R $lineage_ids $gb_summary
    """
}

process download_sequences {
    cpus 8
    memory "64 GB"
    publishDir params.out, mode: 'copy'

    input:
    path(matches)

    output:
    tuple path("sequences/*.fna.gz"), path("sequences/manifest_contigs.csv")

    script:
    """
    download.R --matches $matches \\
      --threads $task.cpus \\
      --out_dir sequences
    
    # Make sure the output file exists, even if empty
    if [ ! -f "sequences/manifest_contigs.csv" ]; then
        echo "Creating empty manifest_contigs.csv file"
        mkdir -p sequences
        echo "id,db,filename,num_records,seqlength" > sequences/manifest_contigs.csv
    fi
    """
}

process download_genome {
    tag "$id"
    cpus 1
    memory "4 GB"
    publishDir "${params.out}/sequences", mode: 'copy'
    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2
    
    input:
    tuple val(id), val(url), val(taxid)
    
    output:
    path "${id}.fna.gz", optional: true
    path "manifest_genome_${id}.csv", optional: true
    
    script:
    """
    echo "Downloading from URL: ${url}"
    
    download_genome.R \\
      --id      ${id} \\
      --url     ${url} \\
      --taxid   ${taxid} \\
      --out_dir .
    
    if [ -f "${id}.fna.gz" ]; then
        echo "Successfully downloaded ${id}"
        # Ensure manifest was created
        if [ ! -f "manifest_genome_${id}.csv" ]; then
            echo "Error: Manifest file not created for ${id}" >&2
            exit 1
        fi
    else
        echo "Warning: Failed to download ${id}" >&2
        exit 1
    fi
    """
}

process make_manifest {
    publishDir "${params.out}/dbs", mode: 'copy'

    input:
    path contigs_meta
    path genome_files

    output:
    path "manifest.csv"

    script:
    """
    Rscript - << 'EOF'
    library(data.table)

    # Load the contigs manifest
    contigs <- fread("${contigs_meta}")

    # Load all genome manifests
    genome_files <- list.files('.', 'manifest_genome_.*\\.csv', full.names=TRUE)
    if (length(genome_files) > 0) {
        genomes <- rbindlist(lapply(genome_files, fread))
    } else {
        genomes <- data.table()
    }

    # Combine and write
    full <- rbind(contigs, genomes, fill=TRUE)
    fwrite(full, 'manifest.csv')
EOF
    """
}

process food_mappings {
    cpus 1
    memory "64 GB"
    publishDir "${params.out}/dbs", mode: 'copy'

    input:
    path(matches)

    output:
    tuple path("food_matches.csv"), path("food_contents.csv.gz")

    script:
    """
    food_mapping.R ${params.out}/dbs/foodb $matches
    """
}

process sketch {
    cpus 2
    memory "4 GB"
    publishDir "${params.out}/sketches", mode: 'copy'

    input:
    path(seq)

    output:
    path("*.sig")

    script:
    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${seq}
    """
}

process ani {
    cpus params.threads
    memory "64 GB"
    publishDir "${params.out}", mode: "copy"

    input:
    path(sigs)

    output:
    path("mash_ani.csv")

    script:
    """
    sourmash compare -k 21 --ani -p ${task.cpus} --csv mash_ani.csv ${sigs}
    """
}