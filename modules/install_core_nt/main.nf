process INSTALL_CORE_NT {

    tag 'core_nt'
    label 'long_serial'
    cache false

    conda "${moduleDir}/environment.yml"

    input:
    val install_dir

    output:
    path 'core_nt.installed', emit: installed

    script:
    """
    set -euo pipefail

    mkdir -p "${install_dir}"

    (
        cd "${install_dir}"

        # Download and decompress core_nt.
        status=0
        update_blastdb.pl --decompress core_nt || status=\$?

        # Exit codes 0 and 1 indicate success.
        if [ "\$status" -gt 1 ]; then
            exit "\$status"
        fi

        # Verify the installed BLAST database.
        blastdbcmd -db core_nt -info > /dev/null
    )

    touch core_nt.installed
    """
}
