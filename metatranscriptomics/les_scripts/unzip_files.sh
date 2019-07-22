
cd $PI_SCRATCH/resilience/metatranscriptomics/raw

for a in *.tgz
do
    a_dir=`expr $a : '\(.*\).tgz'`
    mkdir -p $a_dir
    tar -xvzf $a -C $a_dir
done

