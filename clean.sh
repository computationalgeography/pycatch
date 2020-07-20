rm -f *pyc q.map aaf.map filterRR* piet store1 store2 test.map test.txt multi
rm  filterRR*csv filterSIR*csv *npy
rm -f *npy locats
rm observations/* filterRR*csv filterSIR*csv
for((i = 1; i <= 5000; ++i)); do rm -rf ${i}; done
rm -rf __pycache__
rm -rf pcrasterModules/*.pyc
rm -rf pcrasterModules/__pycache__
