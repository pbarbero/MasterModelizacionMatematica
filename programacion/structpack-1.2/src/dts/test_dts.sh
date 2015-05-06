#! /bin/sh
#./dts --help
export SIZE=16384
export BLOCKSIZE=384
echo "Threads 1, 2, 4, 8, 16 and 32"
echo "No core affinity"
echo "Without pivoting"
OMP_NUM_THREADS=1 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results --raw-header
OMP_NUM_THREADS=2 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=4 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=8 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=16 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=32 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE --raw-results
echo "With pivoting"
OMP_NUM_THREADS=1 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results --raw-header
OMP_NUM_THREADS=2 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=4 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=8 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=16 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results
OMP_NUM_THREADS=32 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE --raw-results
echo " "
echo "Core affinity"
echo "Without pivoting"
OMP_NUM_THREADS=1 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results --raw-header
OMP_NUM_THREADS=2 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=4 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=8 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=16 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=32 OMP_SCHEDULE=static,1 ./dts --size $SIZE --random-seed=1234 -b $BLOCKSIZE -a --raw-results
echo "With pivoting"
OMP_NUM_THREADS=1 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results --raw-header
OMP_NUM_THREADS=2 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=4 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=8 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=16 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results
OMP_NUM_THREADS=32 OMP_SCHEDULE=static,1 ./dts --size $SIZE -p --random-seed=1234 -b $BLOCKSIZE -a --raw-results

STATUS=$?

exit $STATUS
