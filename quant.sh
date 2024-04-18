for sample in `ls data/reads/*_1.fq.gz`
do

dir="data/reads"
base=$(basename $sample "_1.fq.gz")

salmon quant -i salmon_test/concat_index -l A --validateMappings -p 16 -o salmon_test/${base}_quant -1 ${dir}/${base}_1.fq.gz -2 ${dir}/${base}_2.fq.gz
done
