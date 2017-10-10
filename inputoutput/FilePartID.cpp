#include "FilePartID.h"

FilePartID::FilePartID(Collective* col) 
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  const char* filename = col->getTestIDfile().c_str();
  if (rank == 0) printf("ID file: %s\n", filename);

  FILE* fd = fopen(filename, "r");
  if (fd == 0) {
    if (rank == 0) printf("Error: file %s cannot be opened\n", filename);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }

  ns = 0;

  size_t lenmax = 256;
  char* txt;
  char* saveptr;
  char* buffer = (char*) malloc(lenmax*sizeof(char));
  while (getline(&buffer, &lenmax, fd) != -1) {
    if (buffer[0] == '#' || buffer[0] == '\n') continue;
    txt = strtok_r(buffer, " ", &saveptr);
    int s =  atoi(txt);
    txt = strtok_r(NULL, " ", &saveptr);
    long long tR =  atoi(txt);
    txt = strtok_r(NULL, " ", &saveptr);
    long long tP =  atoi(txt);

    bool newspec = true;
    int ind;
    for (int i=0; i<ns; i++) {
      if (s == spec[i]) {
        newspec = false;
        ind = i;
      }
    }

    if (newspec) {
      if (ns == NS_MAX) {
        printf("Too many species in %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
	    ind = ns;
      spec[ind] = s;
      np[ind] = 0;
      tagR[ind] = (long long*) malloc(NP_MAX*sizeof(long long));
      tagP[ind] = (long long*) malloc(NP_MAX*sizeof(long long));
      ns++;
    } 

    tagR[ind][np[ind]] = tR;
    tagP[ind][np[ind]] = tP;
    np[ind]++;
  }

  fclose(fd);
  free(buffer);
}

FilePartID::~FilePartID()
{
  for (int i=0; i<ns; i++) {
     free(tagR[i]);
     free(tagP[i]);
  }
}

void FilePartID::print()
{
  for (int i=0; i<ns; i++) {
    for(int j=0; j<np[i]; j++) {
      printf("%d %ld %ld\n", spec[i], tagR[i][j], tagP[i][j]);
    }
  }
}
