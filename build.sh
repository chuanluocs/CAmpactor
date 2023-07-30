# SamplingCA
cd SamplingCA
make
cd ..

# CAmpactor
sed -i.bak 's/, random_seed      (opt_random_seed)/, random_seed      (1)/g' minisat/core/Solver.cc 
sed -i 's/, phase_saving     (opt_phase_saving)/, phase_saving     (0)/g' minisat/core/Solver.cc 
sed -i.bak 's/PRIi64/ PRIi64 /g' minisat/utils/Options.h
cd minisat
cd core
make Solver.o MROOT=..
cd ..
cd ..
make 
rm ./minisat/core/Solver.cc
mv ./minisat/core/Solver.cc.bak ./minisat/core/Solver.cc
rm ./minisat/utils/Options.h
mv ./minisat/utils/Options.h.bak ./minisat/utils/Options.h
