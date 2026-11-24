process EXTRACT_ANALYSIS_TAXIDS {

    tag "${taxids}"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta), path(accession_taxid)
    path taxdump
    val taxids

    output:
    path "*_*.fasta", emit: fasta, optional: true
    path "analysis_taxids_report.tsv", emit: report
    path "versions.yml", emit: versions

    script:
    """

    # Start the report with a header line
    # Each TaxID adds one row: taxid, name, number of sequences extracted

    echo -e "taxid\\tname\\tsequences" > analysis_taxids_report.tsv

    # The TaxIDs are comma-separated, e.g. 4442,4441
    # tr turns the commas into spaces so the loop gets one TaxID at a time

    for taxid in \$(echo ${taxids} | tr ',' ' '); do

        # Get the scientific name of the TaxID
        # Example: 4442 becomes Camellia_sinensis
        # Replace spaces with underscores for the filename

        name=\$(echo \$taxid | taxonkit lineage -n -L --data-dir ${taxdump} | cut -f2 | tr ' ' '_')

        # Get the taxonomic rank of the TaxID
        # Example: species, genus, family

        rank=\$(echo \$taxid | taxonkit lineage -r -L --data-dir ${taxdump} | cut -f2)

        # If the TaxID is not in the taxdump, the name is empty
        # Write it to the report and go to the next TaxID

        if [ -z "\$name" ]; then
            echo -e "\$taxid\\tnot a valid taxid\\t0" >> analysis_taxids_report.tsv
            continue
        fi

        # Find the TaxID and all TaxIDs below it in the taxonomy tree
        # Example: a genus includes all species belonging to that genus
        # Save all these TaxIDs in taxids.txt

        taxonkit list --ids \$taxid --indent "" --data-dir ${taxdump} > taxids.txt

        # Find accessions associated with the TaxIDs in taxids.txt
        # grep searches for matching TaxIDs in the accession-taxid file
        # cut -f1 keeps only the first column (accession numbers)
        # Save the matching accessions in accessions.txt
        # || true: finding no accessions is not an error, the report will show 0

        grep -wFf taxids.txt ${accession_taxid} | cut -f1 > accessions.txt || true

        # Extract FASTA sequences whose IDs match accessions.txt
        # --id-regexp extracts the accession without the version number
        # Example: AB123456.1 becomes AB123456 for matching
        # Keep the original FASTA headers and complete sequences
        # Save the sequences in a separate FASTA file for each TaxID
        
        seqkit grep --id-regexp '^(\\w+)' -f accessions.txt ${fasta} > \${taxid}_\${name}_\${rank}.fasta

        # Count the sequences extracted (0 if none were found)
        # and add a row to the report

        count=\$(grep -c '>' \${taxid}_\${name}_\${rank}.fasta || true)
        # In the report, write the name with spaces: Camellia_sinensis -> Camellia sinensis
        echo -e "\$taxid\\t\${name//_/ }\\t\$count" >> analysis_taxids_report.tsv

    done


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$(taxonkit version | sed 's/taxonkit v//')
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
