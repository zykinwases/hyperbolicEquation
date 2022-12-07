# module load SpectrumMPI/10.1.0
submit-polus-consecutive:
	g++ -O3 -std=c++11 -fopenmp *.h *.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for i in {1..5} ; do \
			bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 1 $$N 1 out\_1\_$$N\_1.txt ; \
			bsub -n 1 -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 1 $$N 2 out\_1\_$$N\_pi.txt ; \
		done \
	done

submit-polus-parallel-mpi:
	g++ -O3 -std=c++11 -fopenmp *h *.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 2 $$N 1 out\_mpi\_$$p\_$$N\_1.txt ; \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./task3 2 $$N 2 out\_mpi\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done

submit-polus-parallel-omp:
	g++ -O3 -std=c++11 -fopenmp *.h *.cpp -o task3 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm
	for N in 128 256 512 ; do \
		for p in 1 4 8 16 32 ; do \
			for i in {1..5} ; do \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 mpiexec ./task3 3 $$N 1 out\_omp\_$$p\_$$N\_1.txt ; \
				bsub -n $$p -W 00:10 -o /dev/null -e /dev/null OMP_NUM_THREADS=4 mpiexec ./task3 3 $$N 2 out\_omp\_$$p\_$$N\_pi.txt ; \
			done \
		done \
	done