process MAKE_FRAGMENT {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::snapatac2=2.5.3'
    container "biocontainers/snapatac2:2.5.3--py310h4b81fae_0"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fragment.gz"),  emit: fragment
    path  "versions.yml" , emit: versions

    script:
    """
    make_fragment.py $bam ${meta.id}   

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapatac2: \$(python -c "import snapatac2;print(snapatac2.__version__)")
    END_VERSIONS
    """
}