process SNAPATAC2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/snapatac2:2.6.4--py310h28e8315_0"
    containerOptions '--env NUMBA_CACHE_DIR=/tmp'

    input:
    tuple val(meta), path(fragment)
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_gtf), path(gtf)


    output:
    tuple val(meta), path("${meta.id}.fragment.gz"),  emit: fragment
    tuple val(meta), path("*.json"),  emit: json
    tuple val(meta), path("*.png"),  emit: png
    path  "versions.yml" , emit: versions

    script:
    """
    run_snapatac2.py \
     --sample ${meta.id} \
     --fragment $fragment \
     --fasta $fasta \
     --gtf $gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snapatac2: \$(python -c "import snapatac2;print(snapatac2.__version__)")
    END_VERSIONS
    """
}