#!/bin/bash
# Author: Rohini Gadde
# This script profiles the given function in this repo.
# Usage: profile/profile.sh <function>
# Run from the pyslavseq repo!

if [[ $# -ne 1 ]]; then
    echo "Valid functions listed in pyproject.toml. Usage: profile/profile.sh <function>"
    exit 1
fi

root=$(git rev-parse --show-toplevel)
if [[ "${root}" != */pyslavseq ]]; then
    echo "This script must be run from the pyslavseq repository!"
    exit 1
fi

datadir="${root}/profile/data"
outdir="${root}/profile/results"

case "$1" in

    pyslavseq_extract_features)
        BGZ="${datadir}/test.bgz"
        CHROMSIZES="${datadir}/genome.genome"
        FA="${datadir}/genome.fa"

        script="${root}/pyslavseq/features/get_window_features_occupied.py"
        args=(
            --genome_fasta_file $FA --library_3_or_5 3 --occupied --min_mapq 40 --min_ya 20
            --max_yg 15 --chromsizes $CHROMSIZES --window_size 750 --window_step 250
            --min_secondary_mapq 20 $BGZ)
        ;;

    rmdup)
        BAM_IN="${datadir}/no_dedup.bam"
        BAM_OUT="${datadir}/dedup.bam"

        script="${root}/pyslavseq/rmdup/slavseq_rmdup_hts.py"
        args=(--bam $BAM_IN --out $BAM_OUT)
        ;;

    add_tags)
        FA="${datadir}/genome.fa"
        CONSENSUS='ATGTACCCTAAAACTTAGAGTATAATAAA'
        PREFIX_LENGTH=`perl -e 'print length($ENV{CONSENSUS})+2'`
        R1_FLANK_LENGTH=750
        R2_FLANK_LENGTH=$PREFIX_LENGTH
        SOFT_CLIP_LENGTH_THRESHOLD=5

        script="${root}/pyslavseq/tags/add_tags_hts.py"
        args=(
            --genome_fasta_file $FA --prefix_length $PREFIX_LENGTH --consensus $CONSENSUS
            --r1_flank_length $R1_FLANK_LENGTH --r2_flank_length $R2_FLANK_LENGTH
            --soft_clip_length_threshold $SOFT_CLIP_LENGTH_THRESHOLD)
        ;;

    *)
        echo "Valid functions listed in pyproject.toml. Usage: profile/profile.sh <function>"
        exit 1
        ;;
esac

# run cProfile
python -m cProfile -o "${outdir}/tmp_$1.txt" "${script}" "${args[@]}" > /dev/null

echo -e 'sort cumtime\nstats' | python -m pstats "${outdir}/tmp_$1.txt" \
> "${outdir}/$1.txt"

rm ${outdir}/tmp*
