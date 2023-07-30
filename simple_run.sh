echo "c now call SamplingCA to generate an initial PCA"
./SamplingCA/SamplingCA -seed $1 -input_cnf_path $2 -output_testcase_path ./SamplingCA/$3 -t_wise $4
echo "c now call CAmpactor to optimize the initial PCA"
./CAmpactor --seed $1 -i $2 --init_PCA_path SamplingCA/$3 -o $3 --twise $4
rm ./SamplingCA/$3
