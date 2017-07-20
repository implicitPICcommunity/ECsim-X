
#include "Memory.h"

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   ifstream stat_stream("/proc/self/stat",ios_base::in);

   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / (1024.0*1024.0);
   resident_set = rss * page_size_kb / 1024.0;
}

void check_mem_usage()
{

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double local_vm, local_rs;
  double max_vm, max_rs;
  double avg_vm, avg_rs;

  process_mem_usage(local_vm, local_rs);


  MPI_Reduce(&local_vm, &max_vm, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_rs, &max_rs, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_vm, &avg_vm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_rs, &avg_rs, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  avg_vm /= size;
  avg_rs /= size;
  //printf(" memory ussage in proc %d: %.2f / %.2f\n",rank, local_rs, local_vm);

  MPI_Barrier(MPI_COMM_WORLD);
  if (local_rs == max_rs) {
    printf("Memory usage (RES/VIRT): average -> %.1f Mb / %.1f Mb \n"
           "                         maximum -> %.1f Mb / %.1f Mb (rank %d) \n"
            , avg_rs, avg_vm, max_rs, max_vm, rank);
           
  }
}

