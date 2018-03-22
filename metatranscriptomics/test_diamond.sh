BASE_DIR=$SCRATCH/Projects/perturbation_16s
APP_DIR=$SCRATCH/applications/bin/
DIAMOND=$APP_DIR/for_sherlock1/diamond
REF_DIR=$BASE_DIR/data/databases
PYSCRIPT_DIR=$BASE_DIR/metatranscriptomics/pyscripts_edited
INPUT_DIR=$BASE_DIR/data/metatranscriptomics/resilience/input/EBFr_Sub
OUTPUT_DIR=$BASE_DIR/data/metatranscriptomics/resilience/output/EBFr_Sub
input_fwd=M3287_EBFsw_72r_TrM31_1P.fq.gz
input_rev=M3287_EBFsw_72r_TrM31_2P.fq.gz
n_threads=10

base="$(echo $input_fwd | cut -d '_' -f1-3)"
fwd=${input_fwd%.fq.gz}
fwd=${fwd%_001.fastq}
rev=${rev%_001.fastq}
rev=${input_rev%.fq.gz}

start=`date +%s`
$DIAMOND blastx --id 85 --query-cover 65 --min-score 60 \
    --threads $n_threads \
    -d $REF_DIR/nr \
    -q $OUTPUT_DIR/assembled/${base}_contigs_unmapped.fq \
    -o $OUTPUT_DIR/diamond/${base}_nr_contigs.dmdout \
    -f 6 -t $OUTPUT_DIR/dmnd_tmp -k 10

$DIAMOND blastx --id 85 --query-cover 65 --min-score 60 \
    --threads $n_threads \
    -d $REF_DIR/nr \
    -q $OUTPUT_DIR/assembled/${base}_unassembled_unmapped.fq \
    -o $OUTPUT_DIR/diamond/${base}_nr_unassembled.dmdout \
    -f 6 -t $OUTPUT_DIR/dmnd_tmp -k 10

end=`date +%s`
runtime=$(((end-start)/60))
echo "DIAMOND (NR) protein alignment: $runtime min" >> \
    $OUTPUT_DIR/time/${base}_time.log

