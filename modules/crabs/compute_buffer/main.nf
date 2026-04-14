process CRABS_COMPUTE_BUFFER {

    output:
    env(buffersize)

    script:
    """
    buffersize=\$(awk -F'\\t' 'length(\$11) > max { max=length(\$11) } END { print max*2 }' ${params.custom_db})
    """
}