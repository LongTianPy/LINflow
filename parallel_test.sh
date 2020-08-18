source ~/miniconda3/bin/activate
python LINflow.py cleanup linflow_test2
python LINflow.py initiate linflow_test2
# python LINflow.py add_genomes linflow_test2/ -s 1 -i fasta_parallel/ -m fasta_parallel/m1.data;
for M in 1 2 3 4 5
do
python LINflow.py add_genomes linflow_test2/ -s 1 -i fasta_parallel/ -m fasta_parallel/m$M.data &
# sleep 1;
done
