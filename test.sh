#$ -S /bin/bash
for infile in *.krm
do
./knot_representations.py $infile && echo "All good"
done
