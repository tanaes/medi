#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters with remote URLs
params.threads = 20
params.out = "${launchDir}/data"
params.foodb = null       // Default is null, will use remote URL if not specified
params.genbank_summary = null
params.taxdump = null

// Remote URLs (used as fallbacks)
def foodb_remote = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
def genbank_summary_remote = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
def taxdump_remote = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

workflow {
    // Check if provided files exist
    def foodb_exists = params.foodb && file(params.foodb).exists()
    def genbank_summary_exists = params.genbank_summary && file(params.genbank_summary).exists()
    def taxdump_exists = params.taxdump && file(params.taxdump).exists()
    
    // Determine sources (local if exists, otherwise remote)
    def foodb_src = foodb_exists ? params.foodb : foodb_remote
    def genbank_summary_src = genbank_summary_exists ? params.genbank_summary : genbank_summary_remote
    def taxdump_src = taxdump_exists ? params.taxdump : taxdump_remote
    
    // Log the sources being used
    log.info "FooDB source: ${foodb_src}" + (foodb_exists ? " (local file)" : " (remote URL)")
    log.info "GenBank summary source: ${genbank_summary_src}" + (genbank_summary_exists ? " (local file)" : " (remote URL)")
    log.info "Taxdump source: ${taxdump_src}" + (taxdump_exists ? " (local file)" : " (remote URL)")

    //
    // 1) Core DB downloads & matching
    //
    download(foodb_src, genbank_summary_src) | 
        get_taxids
    download_taxa_dbs(taxdump_src)
    get_lineage(get_taxids.out.combine(download_taxa_dbs.out)) | 
        match_taxids

    //
    // 2) Contig (nucleotide) download: one monolithic job
    //
    match_taxids.out | 
        download_sequences

    // Use named output
    contig_seqs = download_sequences.out.fasta_files.flatten()

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

    // Use named output for genomes
    genome_seqs = download_genome.out.fasta_files.collect()

    //
    // 5) Merge contigs + genomes and proceed
    all_seqs = contig_seqs.mix(genome_seqs.flatten())

    // Add a unique identifier to each input file for sketching
    all_seqs_tagged = all_seqs.map { file -> 
        def basename = file.baseName
        return [basename, file]
    }
    
    // Pass the directory where manifests are stored rather than collecting files
    make_manifest(
        download_sequences.out.manifest_file,
        "${params.out}/sequences"  // Directory where genome manifests are published
    )

    sketch(all_seqs_tagged)

    // Collect sketch outputs with unique names for ANI
    ani(sketch.out.collect())

    food_mappings(match_taxids.out)
}


process download {
    cpus 1
    publishDir "${params.out}/dbs", mode: 'copy'

    input:
    val foodb_src
    val genbank_summary_src

    output:
    tuple path("foodb"), path("genbank_summary.tsv")

    script:
    // Use absolute paths for local files
    def foodb_path = foodb_src.startsWith('/') ? foodb_src : 
        (foodb_src.startsWith('http') || foodb_src.startsWith('ftp')) ? foodb_src : "${launchDir}/${foodb_src}"
    def genbank_path = genbank_summary_src.startsWith('/') ? genbank_summary_src : 
        (genbank_summary_src.startsWith('http') || genbank_summary_src.startsWith('ftp')) ? genbank_summary_src : "${launchDir}/${genbank_summary_src}"
    
    """
    # Handle FooDB (local or remote)
    if [[ "${foodb_src}" == http* || "${foodb_src}" == ftp* ]]; then
        echo "Downloading FooDB from: ${foodb_src}"
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${foodb_src} -O foodb.tgz && \\
        tar -xf foodb.tgz && \\
        mv foodb_*_csv foodb
    else
        echo "Using local FooDB file: ${foodb_path}"
        if [[ "${foodb_path}" == *.tar.gz ]]; then
            # It's a tar.gz file
            tar -xf "${foodb_path}"
            mv foodb_*_csv foodb
        elif [[ -d "${foodb_path}" ]]; then
            # It's a directory
            cp -r "${foodb_path}" foodb
        else
            echo "Unsupported FooDB format: ${foodb_path}" >&2
            exit 1
        fi
    fi
    
    # Handle GenBank summary (local or remote)
    if [[ "${genbank_summary_src}" == http* || "${genbank_summary_src}" == ftp* ]]; then
        echo "Downloading GenBank summary from: ${genbank_summary_src}"
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${genbank_summary_src} -O genbank_summary.tsv
    else
        echo "Using local GenBank summary: ${genbank_path}"
        cp "${genbank_path}" genbank_summary.tsv
    fi
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
    val(taxdump_src)

    output:
    path("taxdump")

    script:
    // Use absolute path for file checks
    def taxdump_path = taxdump_src.startsWith('/') ? taxdump_src : "${launchDir}/${taxdump_src}"
    
    """
    # Handle taxdump (local or remote)
    mkdir -p taxdump
    
    # Check if source is a URL or local file
    if [[ "${taxdump_src}" == http* || "${taxdump_src}" == ftp* ]]; then
        # It's a URL
        echo "Downloading taxdump from: ${taxdump_src}"
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \\
            ${taxdump_src} -O taxdump.tar.gz && \\
            tar -xf taxdump.tar.gz --directory taxdump
    else
        # It's a local file - use absolute path
        echo "Using local taxdump: ${taxdump_path}"
        if [[ "${taxdump_path}" == *.tar.gz ]]; then
            # It's a tar.gz file
            tar -xf "${taxdump_path}" --directory taxdump
        elif [[ -d "${taxdump_path}" ]]; then
            # It's a directory
            cp -r "${taxdump_path}"/* taxdump/
        else
            echo "Unsupported taxdump format: ${taxdump_path}" >&2
            exit 1
        fi
    fi
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
    path "sequences/*.fna.gz", emit: fasta_files
    path "sequences/manifest_contigs.csv", emit: manifest_file

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
    path "${id}.fna.gz", emit: fasta_files, optional: true
    path "manifest_genome_${id}.csv", emit: manifest_files, optional: true
    
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
    val genome_manifests_dir

    output:
    path "manifest.csv"

    script:
    """
    make_manifest.R --contigs_manifest "${contigs_meta}" --manifests_dir "${genome_manifests_dir}"
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
    tuple val(id), path(seq)

    output:
    path "${id}.sig"

    script:
    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${seq} -o ${id}.sig
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
