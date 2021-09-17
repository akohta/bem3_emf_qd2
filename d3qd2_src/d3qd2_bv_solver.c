#include "bem3_emf_qd2.h"

int main(int argc,char **argv)
{
  DQD2 qd;

  read_dqd2(argc,argv,&qd);
  print_dqd2(&qd);
  //print_dqd2_mksa(&qd);
  initialize_dqd2(&qd);
  output_node_particles(argv[5],&qd);
  
  solve_bieq(&qd);

  printf("\noutput datafile : %s\n",argv[5]);
  dat_write(argv[5],&qd);

  finalize_dqd2(&qd);

  return 0;
}
