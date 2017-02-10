#PBS -N ej1
#PBS -l nodes=1:ppn=4
#PBS -M js.barbosa10@uniandes.edu.co
#PBS -m abe

module load rocks-openmpi_ib
cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
mpiexec -v -n $NPROCS ./a.out
