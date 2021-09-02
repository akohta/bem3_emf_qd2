/*
 * d3qd2_solve_bieq.c
 *
 *  Created on: Aug 22, 2019
 *      Author: ohta
 */

#include "bem3_emf_qd2.h"

// -- struct def --
typedef struct coef_matrix{
  int type; // coefficient integral type
  int MN; // material number

  // coefficient matrix file name of each domain
  char **tgfn; // G matrix name of each domain
  char **thfn; // H matrix name of each domain

  // derivative coefficient data file name
  char **tdgfn;
  char **tdhfn;
  char **tdffn;

  size_t nn; // total element node number
  // compressed sparse row (CSR) storage file data
  size_t na; // matrix size na*na
  size_t nnz; // number of none zero element
  char *aval; // value
  char *aptr; // pointer
  char *aidx; // index
  // right hand side (rhs) vector
  char *b; //
}CMD;


void solve_bieq(DQD2 *md)
{
  void initalize_cmd(CMD *cm);
  void finalize_cmd(CMD *cm);
  void create_cmatrix(CMD *cm,DQD2 *qd);
  void create_tmatrix_csr(CMD *cm,DQD2 *md);
  void solve_tmatrix_csr(CMD *cm,DQD2 *md);
  void solve_eh_bv(CMD *cm,DQD2 *md);
  void q0_solve_eh_bv(CMD *cm,DQD2 *md);
  void q1_solve_eh_bv(CMD *cm,DQD2 *md);
  void solve_deh_bv(CMD *cm,DQD2 *md);

  time_t start,end,ms,me;
  CMD cm;

  printf("\nsolve modified electromagnetic potential boundary value \n");
  time(&start);

  printf("  coefficient matrix          "); fflush(stdout);
  time(&ms);
  cm.type=PREC_DEF_BE; // precision setting 0:4p GL,1:9p or 7p(triangular) GL, 2:GLN p GL, 3: GHN p GL
  cm.MN=md->MN;
  initalize_cmd(&cm);
  create_cmatrix(&cm,md);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  printf("  solve mEMP boundary value   "); fflush(stdout);
  time(&ms);
  cm.nn=(size_t)md->bd.NN;
  cm.na=cm.nn*8;
  create_tmatrix_csr(&cm,md);
  solve_tmatrix_csr(&cm,md);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  printf("solve electromagnteic field boundary value \n");
  printf("  solve EH boundary value     "); fflush(stdout);
  time(&ms);
  if(SW_DBIEQ==0) q1_solve_eh_bv(&cm,md);
  else solve_eh_bv(&cm,md);
  solve_deh_bv(&cm,md);
  time(&me);
  printf("finished. Elapsed time : %5g (sec)\n",difftime(me,ms));

  finalize_cmd(&cm);
  time(&end);
  printf("Total elapsed time : %g (sec)\n",difftime(end,start));
}

//////////////////////////////////////////////////////////////////////////////////////////
void initalize_cmd(CMD *cm)
{
  int i,MN;

  MN=cm->MN;
  cm->tgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tgfn");
  cm->thfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->thfn");
  cm->tdgfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdgfn");
  cm->tdhfn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdhfn");
  cm->tdffn=(char **)m_alloc2(MN+1,sizeof(char *),"solve_bieq.c, initialize_cmd(), cm->tdffn");
  for(i=0;i<=MN;i++){
    cm->tgfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tgfn[i]");
    cm->thfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->thfn[i]");
    sprintf(cm->tgfn[i],"tmpG_%05d.dat",i);
    sprintf(cm->thfn[i],"tmpH_%05d.dat",i);
    cm->tdgfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdgfn[i]");
    cm->tdhfn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdhfn[i]");
    cm->tdffn[i]=(char *)m_alloc2(16,sizeof(char ),"solve_bieq.c, initialize_cmd(), cm->tdffn[i]");
    sprintf(cm->tdgfn[i],"tmpdG_%05d.dat",i);
    sprintf(cm->tdhfn[i],"tmpdH_%05d.dat",i);
    sprintf(cm->tdffn[i],"tmpdF_%05d.dat",i);
  }

  cm->aval=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aval");
  cm->aptr=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aptr");
  cm->aidx=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->aidx");
  sprintf(cm->aval,"tmpAval.dat");
  sprintf(cm->aptr,"tmpAprt.dat");
  sprintf(cm->aidx,"tmpAidx.dat");
  cm->b=(char *)m_alloc2(16,sizeof(char ),"sovle_bieq.c, initialize_cmd(), cm->b");
  sprintf(cm->b,"tmpB.dat");
}

void finalize_cmd(CMD *cm)
{
  int i;

  // delete temporary file
  for(i=0;i<=cm->MN;i++){
    remove(cm->tgfn[i]);
    remove(cm->thfn[i]);
    remove(cm->tdgfn[i]);
    remove(cm->tdhfn[i]);
    remove(cm->tdffn[i]);
  }
  remove(cm->aval);
  remove(cm->aptr);
  remove(cm->aidx);
  remove(cm->b);

  // free memory
  for(i=0;i<=cm->MN;i++){
    free(cm->tgfn[i]);    free(cm->thfn[i]);
    free(cm->tdgfn[i]);    free(cm->tdhfn[i]);  free(cm->tdffn[i]);
  }
  free(cm->tgfn);  free(cm->thfn);
  free(cm->tdgfn);  free(cm->tdhfn); free(cm->tdffn);

  cm->MN=0;

  free(cm->aval);
  free(cm->aptr);
  free(cm->aidx);
  cm->nn=0;
  cm->na=0;
  cm->nnz=0;
  free(cm->b);
}

void create_cmatrix(CMD *cm,DQD2 *qd)
{
  void create_cmatrix_domain(int did,CMD *cm,DQD2 *qd);

  int i;

  for(i=0;i<=qd->MN;i++) create_cmatrix_domain(i,cm,qd);
}

void create_cmatrix_domain(int did,CMD *cm,DQD2 *qd)
{
  FILE *fg,*fh,*fdg,*fdh,*fdf;
  double complex *tG,*tH,*tdG,*tdH,CC[9],dCz[9],dCe[9],kc;
  double F,vtz[3],vte[3],dFz,dFe,*tdF;
  size_t Ne,N,t,s,tn,tl,i;
  int td,sd;

  Ne=(size_t)qd->bd.sb[did].Ne;
  N=4*Ne;

  tG=(double complex *)m_alloc2(N*N,sizeof(double complex),"create_matrix_domain(),tG"); // malloc
  tH=(double complex *)m_alloc2(N*N,sizeof(double complex),"create_matrix_domain(),tH"); // malloc
  if((fg=fopen(cm->tgfn[did],"wb"))==NULL){    printf("create_matrix_domain(),*fg. Failed to create %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"wb"))==NULL){    printf("create_matrix_domain(),*fh. Failed to create %s file.\n",cm->thfn[did]);    exit(1);  }

  tdG=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"create_matrix_domain(),tdG"); // malloc
  tdH=(double complex *)m_alloc2(2*N*N,sizeof(double complex),"create_matrix_domain(),tdH"); // malloc
  tdF=(double *)m_alloc2(3*N,sizeof(double),"create_matrix_domain(),tdF"); // malloc
  if((fdg=fopen(cm->tdgfn[did],"wb"))==NULL){    printf("create_matrix_domain(),*fdg. Failed to create %s file.\n",cm->tdgfn[did]);    exit(1);  }
  if((fdh=fopen(cm->tdhfn[did],"wb"))==NULL){    printf("create_matrix_domain(),*fdh. Failed to create %s file.\n",cm->tdhfn[did]);    exit(1);  }
  if((fdf=fopen(cm->tdffn[did],"wb"))==NULL){    printf("create_matrix_domain(),*fdf. Failed to create %s file.\n",cm->tdffn[did]);    exit(1);  }


  kc=qd->kn[did];//  printf("did=%d\n",did);
#pragma omp parallel for schedule(dynamic) private(td,tl,tn,vtz,vte,F,dFz,dFe,s,sd,CC,dCz,dCe,i)
  for(t=1;t<=Ne;t++){
    //printf("test t=%zu\n",t);
    td=qd->bd.sb[did].sid[t];
    if( ELT3==check_element_type(td,&(qd->bd)) ) tl=3;
    else tl=4;

    for(tn=0;tn<tl;tn++){
      tz_te_bd_node(vtz,vte,td,tn,&(qd->bd));

      F=0.0;
      dFz=0.0;
      dFe=0.0;
      for(s=1;s<=Ne;s++){
        //printf("test s=%zu\n",s);
        sd=qd->bd.sb[did].sid[s];

        if(did==0 && qd->bd.ps==1){
          if(SW_DBIEQ==0) q0_dcoef_bd_node_t2_cc(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(qd->bd));
          else             q0_dcoef_bd_node_t2   (CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(qd->bd));
        }
        else dcoef_bd_node_t2(CC,dCz,dCe,td,tn,vtz,vte,sd,kc,cm->type,&(qd->bd));

        for(i=0;i<4;i++){
          tG[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+0];
          tH[(t-1)*4*N+tn*N+(s-1)*4+i]=CC[i+4];

          tdG[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*0+(s-1)*4+i]=dCz[i+4];
          tdG[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+0];
          tdH[(t-1)*2*4*N+tn*2*N+N*1+(s-1)*4+i]=dCe[i+4];
        }
        F+=creal(CC[8]);
        dFz+=creal(dCz[8]);
        dFe+=creal(dCe[8]);
      }

      if(did==0){
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=1.0+F;
        tdF[(t-1)*4*3+tn*3+0]=1.0+F;
      }
      else {
        tH[(t-1)*4*N+tn*N+(t-1)*4+tn]+=F;
        tdF[(t-1)*4*3+tn*3+0]=F;
      }
      tdF[(t-1)*4*3+tn*3+1]=dFz;
      tdF[(t-1)*4*3+tn*3+2]=dFe;
    }
  }

  fwrite(tG,sizeof(double complex),N*N,fg);
  fwrite(tH,sizeof(double complex),N*N,fh);
  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);

  fwrite(tdG,sizeof(double complex),2*N*N,fdg);
  fwrite(tdH,sizeof(double complex),2*N*N,fdh);
  fwrite(tdF,sizeof(double),3*N,fdf);
  fclose(fdg);
  fclose(fdh);
  fclose(fdf);
  free(tdG);
  free(tdH);
  free(tdF);
}

void create_tmatrix_csr(CMD *cm,DQD2 *md)
{
  void create_matrix_csr_dac(int did,int cid,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,DQD2 *md);

  FILE *av,*ai,*ap,*b;

  size_t did,cid;

  if((av=fopen(cm->aval,"wb"))==NULL){    printf("create_tmatrix_csr(),*av. Failed to create %s file.\n",cm->aval);    exit(1);  }
  if((ai=fopen(cm->aidx,"wb"))==NULL){    printf("create_tmatrix_csr(),*ai. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((ap=fopen(cm->aptr,"wb"))==NULL){    printf("create_tmatrix_csr(),*ap. Failed to create %s file.\n",cm->aidx);    exit(1);  }
  if((b=fopen(cm->b,"wb"))==NULL){    printf("create_tmatrix_csr(),*b. Failed to create %s file.\n",cm->b);    exit(1);  }

  cm->nnz=0; // initialize nnz
  for(did=0;did<=cm->MN;did++)
    for(cid=0;cid<4;cid++) create_matrix_csr_dac(did, cid,av,ap,ai,b,cm,md);

  fwrite(&(cm->nnz),sizeof(size_t),1,ap);
  //printf("nnz=%zu\n",cm->nnz); // test

  fclose(av);
  fclose(ai);
  fclose(ap);
  fclose(b);
}

void create_matrix_csr_dac(int did,int cid,FILE *av,FILE *ap,FILE *ai,FILE *b,CMD *cm,DQD2 *md)
{

  FILE *fg,*fh;

  double complex *tG,*tH,*tA,tB,k1,k2,delt,beps,sigm;
  double vn[3];
  size_t Ne,t,tn,s,asd,mdid,sdid,etype,l,cc,*ti;
  int td,sd;

  if((fg=fopen(cm->tgfn[did],"rb"))==NULL){    printf("create_tmatrix_csr_dac(),*fg. Failed to open %s file.\n",cm->tgfn[did]);    exit(1);  }
  if((fh=fopen(cm->thfn[did],"rb"))==NULL){    printf("create_tmatrix_csr_dac(),*fh. Failed to open %s file.\n",cm->thfn[did]);    exit(1);  }

  Ne=(size_t)md->bd.sb[did].Ne;
  tG=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"create_tmatrix_csr_dac(),tG"); // malloc
  tH=(double complex *)m_alloc2(Ne*4,sizeof(double complex),"create_tmatrix_csr_dac(),tH"); // malloc

  tA=(double complex *)m_alloc2(cm->na,sizeof(double complex),"create_tmatrix_csr_dac(),tA"); // malloc
  ti=(size_t *)m_alloc2(cm->na,sizeof(size_t),"create_tmatrix_csr_dac(),ti"); // malloc

  for(t=1;t<=Ne;t++){
    td=md->bd.sb[did].sid[t];

    for(tn=0;tn<4;tn++){
      fread(tG,sizeof(double complex),Ne*4,fg);
      fread(tH,sizeof(double complex),Ne*4,fh);
      if( tn==3 && ELT3==check_element_type(td,&(md->bd)) )  continue;

      fwrite(&(cm->nnz),sizeof(size_t),1,ap); // write A pointer
      for(l=0;l<cm->na;l++) tA[l]=0.0;
      tB=0.0;

      for(s=1;s<=Ne;s++){
        sd=md->bd.sb[did].sid[s]; // signed element id
        asd=(size_t)abs(sd);
        mdid=(size_t)md->bd.md[asd]; // main domain id
        sdid=(size_t)md->bd.sd[asd]; // sub domain id
        etype=check_element_type(sd,&(md->bd));        //printf("sd = %d, mdid=%d, sdid=%d\n",sd,mdid,sdid); // test

        if(did==mdid){ // main domain
          if(etype==ELT3){ // linear-triangular element
            for(l=0;l<3;l++) {
              tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
          else { // bi-linear element
            for(l=0;l<4;l++) {
              tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
              tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]=-tG[(s-1)*4+l];
            }
          }
        } // end main domain
        else { // sub domain
          k1=md->kn[mdid];
          k2=md->kn[sdid];
          delt=k1*k1-k2*k2;
          beps=k1*k1/(k2*k2);
          sigm=1.0-beps;

          if(cid<3){ // U
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*(    3) + md->bd.eni[asd][l] ]=delt*vn[cid]*tG[(s-1)*4+l];
                tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]= tG[(s-1)*4+l];
                if(mdid==0) tB+=-tH[(s-1)*4+l]*md->bd.Ui[asd][l][cid]-tG[(s-1)*4+l]*md->bd.dUi[asd][l][cid];
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*(cid+0) + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*(    3) + md->bd.eni[asd][l] ]=delt*vn[cid]*tG[(s-1)*4+l];
                tA[ cm->nn*(cid+4) + md->bd.eni[asd][l] ]= tG[(s-1)*4+l];
                if(mdid==0) tB+=-tH[(s-1)*4+l]*md->bd.Ui[asd][l][cid]-tG[(s-1)*4+l]*md->bd.dUi[asd][l][cid];
              }
            }
          } // end U
          else { // varphi
            if(etype==ELT3){ // linear-triangular element
              for(l=0;l<3;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*0 + md->bd.eni[asd][l] ]= sigm*vn[0]*tG[(s-1)*4+l];
                tA[ cm->nn*1 + md->bd.eni[asd][l] ]= sigm*vn[1]*tG[(s-1)*4+l];
                tA[ cm->nn*2 + md->bd.eni[asd][l] ]= sigm*vn[2]*tG[(s-1)*4+l];
                tA[ cm->nn*3 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*7 + md->bd.eni[asd][l] ]= beps*tG[(s-1)*4+l];
                if(mdid==0) tB+=-sigm*tG[(s-1)*4+l]*(vn[0]*md->bd.Ui[asd][l][0]+vn[1]*md->bd.Ui[asd][l][1]+vn[2]*md->bd.Ui[asd][l][2]);
              }
            }
            else { // bi-linear element
              for(l=0;l<4;l++){
                n_bd_node(vn,sd,l,&(md->bd));
                tA[ cm->nn*0 + md->bd.eni[asd][l] ]= sigm*vn[0]*tG[(s-1)*4+l];
                tA[ cm->nn*1 + md->bd.eni[asd][l] ]= sigm*vn[1]*tG[(s-1)*4+l];
                tA[ cm->nn*2 + md->bd.eni[asd][l] ]= sigm*vn[2]*tG[(s-1)*4+l];
                tA[ cm->nn*3 + md->bd.eni[asd][l] ]= tH[(s-1)*4+l];
                tA[ cm->nn*7 + md->bd.eni[asd][l] ]= beps*tG[(s-1)*4+l];
                if(mdid==0) tB+=-sigm*tG[(s-1)*4+l]*(vn[0]*md->bd.Ui[asd][l][0]+vn[1]*md->bd.Ui[asd][l][1]+vn[2]*md->bd.Ui[asd][l][2]);
              }
            }
          } // end varphi
        } // end sub domain
      } // end for s

      // compress and store data
      cc=0;
      for(l=0;l<cm->na;l++){
        if( creal(tA[l])==0.0 && cimag(tA[l])==0.0) continue;
        tA[cc]=tA[l];
        ti[cc]=l;
        cc+=1;
      }
      fwrite(tA,sizeof(double complex),cc,av);
      fwrite(ti,sizeof(size_t),cc,ai);
      cm->nnz+=cc;

      fwrite(&tB,sizeof(double complex),1,b);

    } // end for tn

  } // end for t

  fclose(fg);
  fclose(fh);
  free(tG);
  free(tH);
  free(tA);
  free(ti);
}

void solve_tmatrix_csr(CMD *cm,DQD2 *md)
{
  int b_mkl_solver_pardiso_d(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solver.c
  int b_mkl_solver_pardiso_s(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solver.c
  void tmatrix_bd_store(double complex *X,DQD2 *md);

  int err;

  double complex *x; // results
  x=(double complex *)m_alloc2(cm->na,sizeof(double complex),"solve_tmatrix_csr(),x");

  if(PARDISO_PREC==0) err=b_mkl_solver_pardiso_d(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  else err=b_mkl_solver_pardiso_s(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  if(err!=0){
    printf("b_mkl_solver_pardiso() returned error. err=%d. Exit...\n",err);
    exit(1);
  }

  /*
  // for test
  int b_mkl_solver_lapacke(size_t n,size_t nnz,char *fn_a,char *fn_xa,char *fn_asub,char *fn_b,double complex *x); // b_mkl_solfer.c
  err=b_mkl_solver_lapacke(cm->na,cm->nnz,cm->aval,cm->aptr,cm->aidx,cm->b,x);
  for(i=10;i<15;i++){
    printf("x[%d]=%15.14e + %15.14e\n",i,creal(x[i]),cimag(x[i]));
  }
 */

  // store data
  tmatrix_bd_store(x,md);

  free(x);
}

void tmatrix_bd_store(double complex *X,DQD2 *md)
{
  double complex k1,k2,delt,beps,sigm,tU[4],tUi[4];
  double n[3];
  size_t d,s,l,nn,i;
  int sd,asd,etype,mdid,sdid;

  nn=md->bd.NN;

  for(d=0;d<=md->MN;d++){
    for(s=1;s<=md->bd.sb[d].Ne;s++){
      sd=md->bd.sb[d].sid[s];
      asd=abs(sd);
      mdid=md->bd.md[asd];
      sdid=md->bd.sd[asd];
      etype=check_element_type(sd,&(md->bd));

      if(d==mdid){ // main domain
        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            for(i=0;i<4;i++){ // store U, dU
              md->bd.sb[d]. U[s][l][i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d].dU[s][l][i]=X[nn*(i+4)+md->bd.eni[asd][l]];
            }
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            for(i=0;i<4;i++){ // store U, dU
              md->bd.sb[d]. U[s][l][i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d].dU[s][l][i]=X[nn*(i+4)+md->bd.eni[asd][l]];
            }
          }
        }
      }
      else { // subdomain
        k1=md->kn[mdid];
        k2=md->kn[sdid];
        delt=k1*k1-k2*k2;
        beps=k1*k1/(k2*k2);
        sigm=1.0-beps;

        if(etype==ELT3){ // linear triangular element
          for(l=0;l<3;l++){ // node id
            n_bd_node(n,asd,l,&(md->bd));
            tU[3]=X[nn*3+md->bd.eni[asd][l]];
            for(i=0;i<3;i++){ // store U, dU
              tU[i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d]. U[s][l][i]=tU[i];
              md->bd.sb[d].dU[s][l][i]=-(X[nn*(i+4)+md->bd.eni[asd][l]]+delt*tU[3]*n[i]);
              if(mdid==0){
                tUi[i]=md->bd.Ui[asd][l][i];
                md->bd.sb[d]. U[s][l][i]+= tUi[i];
                md->bd.sb[d].dU[s][l][i]+=-md->bd.dUi[asd][l][i];
              }
            }
            // store phi, dphi
            md->bd.sb[d]. U[s][l][3]=tU[3];
            md->bd.sb[d].dU[s][l][3]=-(beps*X[nn*7+md->bd.eni[asd][l]]+sigm*(n[0]*tU[0]+n[1]*tU[1]+n[2]*tU[2]));
            if(mdid==0)  md->bd.sb[d].dU[s][l][3]+=-sigm*(n[0]*tUi[0]+n[1]*tUi[1]+n[2]*tUi[2]);
          }
        }
        else { // bi-linear element
          for(l=0;l<4;l++){ // node id
            n_bd_node(n,asd,l,&(md->bd));
            tU[3]=X[nn*3+md->bd.eni[asd][l]];
            for(i=0;i<3;i++){ // store U, dU
              tU[i]=X[nn*(i+0)+md->bd.eni[asd][l]];
              md->bd.sb[d]. U[s][l][i]=tU[i];
              md->bd.sb[d].dU[s][l][i]=-(X[nn*(i+4)+md->bd.eni[asd][l]]+delt*tU[3]*n[i]);
              if(mdid==0){
                tUi[i]=md->bd.Ui[asd][l][i];
                md->bd.sb[d]. U[s][l][i]+= tUi[i];
                md->bd.sb[d].dU[s][l][i]+=-md->bd.dUi[asd][l][i];
              }
            }
            // store phi, dphi
            md->bd.sb[d]. U[s][l][3]=tU[3];
            md->bd.sb[d].dU[s][l][3]=-(beps*X[nn*7+md->bd.eni[asd][l]]+sigm*(n[0]*tU[0]+n[1]*tU[1]+n[2]*tU[2]));
            if(mdid==0)  md->bd.sb[d].dU[s][l][3]+=-sigm*(n[0]*tUi[0]+n[1]*tUi[1]+n[2]*tUi[2]);
          }
        }
      } // end subdomain

    }  // end s

  } // end d

}

void solve_eh_bv(CMD *cm,DQD2 *md)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dUz[4],dUe[4],ch;
  double tdF[3],vtz[3],vte[3],vn[3],a[9],sig;
  size_t Ne,d,t,atd,tn,i;
  int at;

  ch=md->pw.lambda_0/(2.0*M_PI*I);

  for(d=0;d<=md->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=md->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_eh_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_eh_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=md->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(md->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*md->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dUdtz,dUdte
        dudt_fcm(dUz,dUe,t,tn,Ne,md->bd.sb[d].U,md->bd.sb[d].dU,tdG,tdH,tdF);
        // electric field
        for(i=0;i<3;i++) md->bd.sb[d].E[t][tn][i]=-(md->bd.sb[d].dU[t][tn][3]*a[i*3+0]+dUz[3]*a[i*3+1]+dUe[3]*a[i*3+2])+md->bd.sb[d].U[t][tn][i];
/*
        // test E
        printf("d=%zu,t=%zu,tn=%zu\n",d,t,tn);
        for(i=0;i<3;i++) printf("E=%15.14e + %15.14e\n",creal(md->bd.sb[d].E[t][tn][i]),cimag(md->bd.sb[d].E[t][tn][i]));
        double complex t0,t1;
        t0=vn[0]*md->bd.sb[d].E[t][tn][0]+vn[1]*md->bd.sb[d].E[t][tn][1]+vn[2]*md->bd.sb[d].E[t][tn][2];
        t1=-md->bd.sb[d].dU[t][tn][3]+vn[0]*md->bd.sb[d].U[t][tn][0]+vn[1]*md->bd.sb[d].U[t][tn][1]+vn[2]*md->bd.sb[d].U[t][tn][2];
        printf("t0=%15.14e + %15.14e\n",creal(t0),cimag(t0));
        printf("t1=%15.14e + %15.14e\n",creal(t1),cimag(t1));
*/
        // magnetic field
        md->bd.sb[d].H[t][tn][0]=ch*( (md->bd.sb[d].dU[t][tn][2]*a[1*3+0]+dUz[2]*a[1*3+1]+dUe[2]*a[1*3+2])
                                     -(md->bd.sb[d].dU[t][tn][1]*a[2*3+0]+dUz[1]*a[2*3+1]+dUe[1]*a[2*3+2]) );
        md->bd.sb[d].H[t][tn][1]=ch*( (md->bd.sb[d].dU[t][tn][0]*a[2*3+0]+dUz[0]*a[2*3+1]+dUe[0]*a[2*3+2])
                                     -(md->bd.sb[d].dU[t][tn][2]*a[0*3+0]+dUz[2]*a[0*3+1]+dUe[2]*a[0*3+2]) );
        md->bd.sb[d].H[t][tn][2]=ch*( (md->bd.sb[d].dU[t][tn][1]*a[0*3+0]+dUz[1]*a[0*3+1]+dUe[1]*a[0*3+2])
                                     -(md->bd.sb[d].dU[t][tn][0]*a[1*3+0]+dUz[0]*a[1*3+1]+dUe[0]*a[1*3+2]) );
      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d
}

void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF)
{
  size_t i,s,sn;

  for(i=0;i<4;i++){
    dUz[i]=0.0;
    dUe[i]=0.0;
  }
  for(s=1;s<=Ne;s++){
    for(sn=0;sn<4;sn++){
      for(i=0;i<4;i++){
        dUz[i]+=tdG[Ne*4*0+4*(s-1)+sn]*dU[s][sn][i]-tdH[Ne*4*0+4*(s-1)+sn]*U[s][sn][i];
        dUe[i]+=tdG[Ne*4*1+4*(s-1)+sn]*dU[s][sn][i]-tdH[Ne*4*1+4*(s-1)+sn]*U[s][sn][i];
      }
    }
  }
  for(i=0;i<4;i++){
    dUz[i]=(dUz[i]-U[t][tn][i]*tdF[1])/tdF[0];
    dUe[i]=(dUe[i]-U[t][tn][i]*tdF[2])/tdF[0];
  }
}

void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte)
{
  double i_det;
  int i;

  i_det=1.0/(vn[0]*vtz[1]*vte[2]+vn[1]*vtz[2]*vte[0]+vn[2]*vtz[0]*vte[1]
          -( vn[2]*vtz[1]*vte[0]+vn[1]*vtz[0]*vte[2]+vn[0]*vtz[2]*vte[1]));

  ma[0*3+0]= vtz[1]*vte[2]-vtz[2]*vte[1];
  ma[0*3+1]=- vn[1]*vte[2]+ vn[2]*vte[1];
  ma[0*3+2]=  vn[1]*vtz[2]- vn[2]*vtz[1];

  ma[1*3+0]=-vtz[0]*vte[2]+vtz[2]*vte[0];
  ma[1*3+1]=  vn[0]*vte[2]- vn[2]*vte[0];
  ma[1*3+2]=- vn[0]*vtz[2]+ vn[2]*vtz[0];

  ma[2*3+0]= vtz[0]*vte[1]-vtz[1]*vte[0];
  ma[2*3+1]=- vn[0]*vte[1]+ vn[1]*vte[0];
  ma[2*3+2]=  vn[0]*vtz[1]- vn[1]*vtz[0];

  for(i=0;i<9;i++) ma[i]*=i_det;
}

void solve_deh_bv(CMD *cm,DQD2 *md)
{
  FILE *fg,*fh;
  MKL_Complex16 *A,*B;
  double complex *tG,*tH;
  size_t Ne,d,t,atd,tn,nn,nc,nr,s,ns,asd,i;
  MKL_INT *ipiv,nrhs,lda,ldb,info;

  nrhs=6;

  for(d=0;d<=md->MN;d++){
    Ne=md->bd.sb[d].Ne;
    if(Ne==0) continue;

    if((fg=fopen(cm->tgfn[d],"rb"))==NULL){    printf("solve_deh_bv(),*fg. Failed to open %s file.\n",cm->tgfn[d]);    exit(1);  }
    if((fh=fopen(cm->thfn[d],"rb"))==NULL){    printf("solve_deh_bv(),*fh. Failed to open %s file.\n",cm->thfn[d]);    exit(1);  }

    // matrix size
    nn=0;
    for(s=1;s<=Ne;s++){
      asd=abs(md->bd.sb[d].sid[s]);
      if( ELT3==check_element_type(asd,&(md->bd)) ) nn+=3;
      else nn+=4;
    }
    //printf("nn=%zu\n",nn); // test
    lda=nn;
    ldb=6;

    A =(MKL_Complex16 *)m_alloc2(nn*nn,sizeof(MKL_Complex16),"solve_deh_bv(),A"); // malloc
    B =(MKL_Complex16 *)m_alloc2(nn*6,sizeof(MKL_Complex16),"solve_dev_bv(),B"); // malloc
    ipiv=(MKL_INT *)m_alloc2(nn,sizeof(MKL_INT),"solve_deh_bv(),ipiv"); // malloc
    tG=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"solve_dev_bv(),tG"); // malloc
    tH=(double complex *)m_alloc2(4*Ne,sizeof(double complex),"solve_dev_bv(),tH"); // malloc

    nc=0;
    for(t=1;t<=Ne;t++){
      atd=abs(md->bd.sb[d].sid[t]);

      for(tn=0;tn<4;tn++){
        fread(tG,sizeof(double complex),Ne*4,fg);
        fread(tH,sizeof(double complex),Ne*4,fh);
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        for(i=0;i<6;i++){
          B[6*nc+i].real=0.0;
          B[6*nc+i].imag=0.0;
        }
        nr=0;
        for(s=1;s<=Ne;s++){
          asd=abs(md->bd.sb[d].sid[s]);
          if( ELT4==check_element_type(asd,&(md->bd)) ){
            for(ns=0;ns<4;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[6*nc+i+0].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+0].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+3].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
                B[6*nc+i+3].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
              }
            }
          }
          else {
            for(ns=0;ns<3;ns++) {
              A[nn*nc+nr].real=creal(tG[4*(s-1)+ns]);
              A[nn*nc+nr].imag=cimag(tG[4*(s-1)+ns]);
              nr+=1;
              for(i=0;i<3;i++){
                B[6*nc+i+0].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+0].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].E[s][ns][i]);
                B[6*nc+i+3].real+=creal(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
                B[6*nc+i+3].imag+=cimag(tH[4*(s-1)+ns]*md->bd.sb[d].H[s][ns][i]);
              }
            }
          }
        } // end for s
        //printf("nr=%zu\n",nr);
        nc+=1;
      } // end for tn

    } // end for t
    //printf("nc=%zu\n",nc);

    // solve mkl lapack
    info = LAPACKE_zgesv( LAPACK_ROW_MAJOR, nn, nrhs, A, lda, ipiv, B, ldb );
    if(info!=0){
      printf("solve_deh_bv(), LAPACKE_zgesv(). info=%lld error. Exit...\n",info);
      exit(1);
    }

    // store data
    nc=0;
    for(s=1;s<=Ne;s++){
      asd=abs(md->bd.sb[d].sid[s]);
      if( ELT4==check_element_type(asd,&(md->bd)) ){
        for(ns=0;ns<4;ns++){
          for(i=0;i<3;i++) {
            md->bd.sb[d].dE[s][ns][i]=B[6*nc+i+0].real+B[6*nc+i+0].imag*I;
            md->bd.sb[d].dH[s][ns][i]=B[6*nc+i+3].real+B[6*nc+i+3].imag*I;
          }
          nc+=1;
        }
      }
      else {
        for(ns=0;ns<3;ns++){
          for(i=0;i<3;i++) {
            md->bd.sb[d].dE[s][ns][i]=B[6*nc+i+0].real+B[6*nc+i+0].imag*I;
            md->bd.sb[d].dH[s][ns][i]=B[6*nc+i+3].real+B[6*nc+i+3].imag*I;
          }
          nc+=1;
        }
      }
    }

    fclose(fg);
    fclose(fh);
    free(A);
    free(B);
    free(ipiv);
    free(tG);
    free(tH);
  } // end for d
}

void q0_solve_eh_bv(CMD *cm,DQD2 *md)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void inv_matrix_c33(double complex *ma,double complex *vn,double complex *vtz,double complex *vte);
  void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dUz[4],dUe[4],ch,*tE0,*tH0,es,ew,Ei[3],Hi[3],dE,rhs[3],esvn[3],ac[9],vtzc[3],vtec[3];
  double tdF[3],vtz[3],vte[3],vn[3],a[9],sig;
  size_t Ne,d,t,atd,tn,i;
  int at;

  ch=md->pw.lambda_0/(2.0*M_PI*I);
  tE0=(double complex *)m_alloc2(md->bd.NN*3,sizeof(double complex),"q0_solve_eh_bv(),tE0"); // malloc
  tH0=(double complex *)m_alloc2(md->bd.NN*3,sizeof(double complex),"q0_solve_eh_bv(),tH0"); // malloc
  es=md->n[0]*md->n[0];

  for(d=1;d<=md->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("q0_solve_eh_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("q0_solve_eh_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("q0_solve_eh_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=md->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"q0_solve_eh_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"q0_solve_eh_bv(),tdH"); // malloc
    ew=md->n[d]*md->n[d];

    for(t=1;t<=Ne;t++){
      at=md->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(md->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*md->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dUdtz,dUdte
        dudt_fcm(dUz,dUe,t,tn,Ne,md->bd.sb[d].U,md->bd.sb[d].dU,tdG,tdH,tdF);
        // electric field
        for(i=0;i<3;i++) md->bd.sb[d].E[t][tn][i]=-(md->bd.sb[d].dU[t][tn][3]*a[i*3+0]+dUz[3]*a[i*3+1]+dUe[3]*a[i*3+2])+md->bd.sb[d].U[t][tn][i];

        // magnetic field
        md->bd.sb[d].H[t][tn][0]=ch*( (md->bd.sb[d].dU[t][tn][2]*a[1*3+0]+dUz[2]*a[1*3+1]+dUe[2]*a[1*3+2])
                                     -(md->bd.sb[d].dU[t][tn][1]*a[2*3+0]+dUz[1]*a[2*3+1]+dUe[1]*a[2*3+2]) );
        md->bd.sb[d].H[t][tn][1]=ch*( (md->bd.sb[d].dU[t][tn][0]*a[2*3+0]+dUz[0]*a[2*3+1]+dUe[0]*a[2*3+2])
                                     -(md->bd.sb[d].dU[t][tn][2]*a[0*3+0]+dUz[2]*a[0*3+1]+dUe[2]*a[0*3+2]) );
        md->bd.sb[d].H[t][tn][2]=ch*( (md->bd.sb[d].dU[t][tn][1]*a[0*3+0]+dUz[1]*a[0*3+1]+dUe[1]*a[0*3+2])
                                     -(md->bd.sb[d].dU[t][tn][0]*a[1*3+0]+dUz[0]*a[1*3+1]+dUe[0]*a[1*3+2]) );

        // for scattered field
        if(md->bd.md[atd]==0){
          for(i=0;i<3;i++) rhs[i]=0.0;
          calc_mipw_EH(Ei,Hi,md->bd.ren[atd][tn],&(md->pw));

          for(i=0;i<3;i++){
            esvn[i]=es*vn[i];
            vtzc[i]=vtz[i];
            vtec[i]=vte[i];
            dE=md->bd.sb[d].E[t][tn][i]-Ei[i];
            rhs[0]+=(ew*md->bd.sb[d].E[t][tn][i]-es*Ei[i])*vn[i];
            rhs[1]+=dE*vtz[i];
            rhs[2]+=dE*vte[i];
          }
          inv_matrix_c33(ac,esvn,vtzc,vtec);

          for(i=0;i<3;i++){
            tE0[md->bd.eni[atd][tn]*3+i]=ac[3*i+0]*rhs[0]+ac[3*i+1]*rhs[1]+ac[3*i+2]*rhs[2];
            tH0[md->bd.eni[atd][tn]*3+i]=md->bd.sb[d].H[t][tn][i]-Hi[i];
          }
        }

      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d

  // store scattered data
  for(t=1;t<=md->bd.sb[0].Ne;t++){
    at=md->bd.sb[0].sid[t];
    atd=abs(at);
    for(tn=0;tn<4;tn++){
      if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;
      for(i=0;i<3;i++){
        md->bd.sb[0].E[t][tn][i]=tE0[md->bd.eni[atd][tn]*3+i];
        md->bd.sb[0].H[t][tn][i]=tH0[md->bd.eni[atd][tn]*3+i];
      }
    }
  }

  free(tE0);
  free(tH0);
}

void inv_matrix_c33(double complex *ma,double complex *vn,double complex *vtz,double complex *vte)
{
  double complex i_det;
  int i;

  i_det=1.0/(vn[0]*vtz[1]*vte[2]+vn[1]*vtz[2]*vte[0]+vn[2]*vtz[0]*vte[1]
          -( vn[2]*vtz[1]*vte[0]+vn[1]*vtz[0]*vte[2]+vn[0]*vtz[2]*vte[1]));

  ma[0*3+0]= vtz[1]*vte[2]-vtz[2]*vte[1];
  ma[0*3+1]=- vn[1]*vte[2]+ vn[2]*vte[1];
  ma[0*3+2]=  vn[1]*vtz[2]- vn[2]*vtz[1];

  ma[1*3+0]=-vtz[0]*vte[2]+vtz[2]*vte[0];
  ma[1*3+1]=  vn[0]*vte[2]- vn[2]*vte[0];
  ma[1*3+2]=- vn[0]*vtz[2]+ vn[2]*vtz[0];

  ma[2*3+0]= vtz[0]*vte[1]-vtz[1]*vte[0];
  ma[2*3+1]=- vn[0]*vte[1]+ vn[1]*vte[0];
  ma[2*3+2]=  vn[0]*vtz[1]- vn[1]*vtz[0];

  for(i=0;i<9;i++) ma[i]*=i_det;
}

void q1_solve_eh_bv(CMD *cm,DQD2 *md)
{
  void inv_matrix_33(double *ma,double *vn,double *vtz,double *vte);
  void dudt_fcm(double complex *dUz,double complex *dUe,int t,int tn,size_t Ne,double complex ***U,double complex ***dU,double complex *tdG,double complex *tdH,double *tdF);

  FILE *fdg,*fdh,*fdf;
  double complex *tdG,*tdH,dUz[4],dUe[4],ch,dUiz[3],dUie[3],e[3],h[3],dh[3];
  double tdF[3],vtz[3],vte[3],vn[3],a[9],sig;
  size_t Ne,d,t,atd,tn,s,i;
  int at;

  ch=md->pw.lambda_0/(2.0*M_PI*I);

  for(d=1;d<=md->MN;d++){
    if((fdg=fopen(cm->tdgfn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdg. Failed to open %s file.\n",cm->tdgfn[d]);    exit(1);  }
    if((fdh=fopen(cm->tdhfn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdh. Failed to open %s file.\n",cm->tdhfn[d]);    exit(1);  }
    if((fdf=fopen(cm->tdffn[d],"rb"))==NULL){    printf("solve_eh_bv(),*fdf. Failed to open %s file.\n",cm->tdffn[d]);    exit(1);  }

    Ne=md->bd.sb[d].Ne;
    tdG=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_eh_bv(),tdG"); // malloc
    tdH=(double complex *)m_alloc2(2*Ne*4,sizeof(double complex),"solve_eh_bv(),tdH"); // malloc

    for(t=1;t<=Ne;t++){
      at=md->bd.sb[d].sid[t];
      if(at<0) sig=-1.0;
      else sig=1.0;
      atd=abs(at);

      for(tn=0;tn<4;tn++){
        fread(tdG,sizeof(double complex),2*Ne*4,fdg);
        fread(tdH,sizeof(double complex),2*Ne*4,fdh);
        fread(tdF,sizeof(double),3,fdf);
        if( tn==3 && ELT3==check_element_type(atd,&(md->bd)) )  continue;

        // tangential vector
        tz_te_bd_node(vtz,vte,atd,tn,&(md->bd));
        // normal vector
        for(i=0;i<3;i++) vn[i]=sig*md->bd.wen[atd][tn][i];
        vuni_d(vn);
        // inv-matrix
        inv_matrix_33(a,vn,vtz,vte);
        // dUdtz,dUdte
        dudt_fcm(dUz,dUe,t,tn,Ne,md->bd.sb[d].U,md->bd.sb[d].dU,tdG,tdH,tdF);
        // electric field
        for(i=0;i<3;i++) md->bd.sb[d].E[t][tn][i]=-(md->bd.sb[d].dU[t][tn][3]*a[i*3+0]+dUz[3]*a[i*3+1]+dUe[3]*a[i*3+2])
                                                    +md->bd.sb[d].U[t][tn][i];
        // magnetic field
        md->bd.sb[d].H[t][tn][0]=ch*( (md->bd.sb[d].dU[t][tn][2]*a[1*3+0]+dUz[2]*a[1*3+1]+dUe[2]*a[1*3+2])
                                     -(md->bd.sb[d].dU[t][tn][1]*a[2*3+0]+dUz[1]*a[2*3+1]+dUe[1]*a[2*3+2]) );
        md->bd.sb[d].H[t][tn][1]=ch*( (md->bd.sb[d].dU[t][tn][0]*a[2*3+0]+dUz[0]*a[2*3+1]+dUe[0]*a[2*3+2])
                                     -(md->bd.sb[d].dU[t][tn][2]*a[0*3+0]+dUz[2]*a[0*3+1]+dUe[2]*a[0*3+2]) );
        md->bd.sb[d].H[t][tn][2]=ch*( (md->bd.sb[d].dU[t][tn][1]*a[0*3+0]+dUz[1]*a[0*3+1]+dUe[1]*a[0*3+2])
                                     -(md->bd.sb[d].dU[t][tn][0]*a[1*3+0]+dUz[0]*a[1*3+1]+dUe[0]*a[1*3+2]) );

        // for scattered field
        if(md->bd.md[atd]==0){
          calc_mipw_EH_dEHdv(e,h,dUiz,dh,md->bd.ren[atd][tn],vtz,&(md->pw));
          calc_mipw_EH_dEHdv(e,h,dUie,dh,md->bd.ren[atd][tn],vte,&(md->pw));

          // scan open region id
          for(s=1;s<=md->bd.sb[0].Ne;s++){
            if(atd==md->bd.sb[0].sid[s]) break;
          }
          if(s>md->bd.sb[0].Ne){
            printf("solve_bidq.c, q1_solve_eh_bv(),open region id error! id=%zu. Exit...\n",s);
            exit(0);
          }

          // electric field
          for(i=0;i<3;i++) md->bd.sb[0].E[s][tn][i]=-(-md->bd.sb[0].dU[s][tn][3]*a[i*3+0]+dUz[3]*a[i*3+1]+dUe[3]*a[i*3+2])
                                                      +md->bd.sb[0].U[s][tn][i];

          // magnetic field
          md->bd.sb[0].H[s][tn][0]=ch*( (-md->bd.sb[0].dU[s][tn][2]*a[1*3+0]+(dUz[2]-dUiz[2])*a[1*3+1]+(dUe[2]-dUie[2])*a[1*3+2])
                                       -(-md->bd.sb[0].dU[s][tn][1]*a[2*3+0]+(dUz[1]-dUiz[1])*a[2*3+1]+(dUe[1]-dUie[1])*a[2*3+2]) );
          md->bd.sb[0].H[s][tn][1]=ch*( (-md->bd.sb[0].dU[s][tn][0]*a[2*3+0]+(dUz[0]-dUiz[0])*a[2*3+1]+(dUe[0]-dUie[0])*a[2*3+2])
                                       -(-md->bd.sb[0].dU[s][tn][2]*a[0*3+0]+(dUz[2]-dUiz[2])*a[0*3+1]+(dUe[2]-dUie[2])*a[0*3+2]) );
          md->bd.sb[0].H[s][tn][2]=ch*( (-md->bd.sb[0].dU[s][tn][1]*a[0*3+0]+(dUz[1]-dUiz[1])*a[0*3+1]+(dUe[1]-dUie[1])*a[0*3+2])
                                       -(-md->bd.sb[0].dU[s][tn][0]*a[1*3+0]+(dUz[0]-dUiz[0])*a[1*3+1]+(dUe[0]-dUie[0])*a[1*3+2]) );

        } // end scattered field

      } // end for tn
    } // end for t
    fclose(fdg);
    fclose(fdh);
    fclose(fdf);
    free(tdG);
    free(tdH);
  } // end for d
}
