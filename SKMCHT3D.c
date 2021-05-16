/*******************************************************************************
Subvolume KMC Simulations of freely diffusing bacteria in a 3D porous media which
causes it to exhibit a two state motion of hopping and trapping.

USAGE

%cc -O3 SKMCHT3D.c -o SKMCHT3D
%./SKMCHT3D > SKMCHT3D.out
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "minHeap.h"

#define SNAPSHOT_RATE 1 /* we take a snapshot after every SNAPSHOT_RATE events */
#define NUM_CELLS_ONESIDE 50 /* number of lattice cells along one side of square system domain */
#define TOTAL_LATTICE_CELLS (NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE) /* Total number of lattice cells */
#define TOTAL_BACTERIA_NUM 200 /* Total number of emerin monomer proteins */
#define N0 (TOTAL_BACTERIA_NUM/TOTAL_LATTICE_CELLS) /* Initial number of proteins in each lattice */
#define V (6.0*6.0*6.0) /* Total area of system domain (LATTICE_CELL_LENGTH^2) */
#define LATTICE_CELL_LENGTH 1.0 /* lattice cell length (=1 micrometer/SCALE)*/
#define TAU_TRAPPED (0.5/2.0)/* time between hops inside nanodomain (s) */
#define TAU_HOP 0.0 /* time between hops outside nanodomain (s). we set to zero since none were reported in the experiments. */
#define SCALE 1 /*re-scales the size of the lattice cells */
#define HOP_AVG (3.24*SCALE) /*average hop length */
#define TIME_MAX 20.0 /*simulation time limit */
#define ITER_MAX 1 /*number of iterations of the simulation */
#define DIM 3 /* Our domain is 3D */

int main() {
  /* directionID will hold lattice neighbor index protein hops into */
  /* lambda will hold the index of the lattice site index it hops out of */
  /* gamma will hold the index of the lattice site index it hops into */
  /* hop will hold sampled hop length */
  /* x,y,z will hold coordinates of the bacteria as it travels from trapped state to trapped state (for 1 bacteria simulation) */
  /* u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11 are to hold randoms numbers */
  /* tt will hold sampled trapped time */
  /* th will hold sampled hop time */
  /* D will hold the measured simulated diffusivity */
  int directionID, lambda, gamma, u4, hop, y, x, z;
  double u1, u2, u3, u5, u6, u7, u8, u9, u10, u11, tt, th, D;

  /* Ntot is for the total number of proteins across 
     all lattice cells */
  int Ntot = 0;

  /* tlambda will hold the sampled event time in lattice lambda */
  /* tgamma will hold the sampled event time in lattice gamma */
  /* Dsum will hold the diffusivities summed over iterations for averaging in later step */
  /* Davg will hold the diffusivity averaged over iterations */
  double tlambda = 0.0, tgamma = 0.0, Dsum = 0.0, Davg = 0.0;

  int N[TOTAL_LATTICE_CELLS]; /* number of emerin monomer proteins in each lattice cell */
  int nn[TOTAL_LATTICE_CELLS][6]; /* nearest neighbors of each lattice cell, we use 
                   periodic boundary conditions at system boundary */
  double W[TOTAL_LATTICE_CELLS]; /* transition rates (diffusion intensities) in each lattice cell */
  double t[TOTAL_LATTICE_CELLS]; /* next queued (sampled) event time in each lattice cell */
  double tau[TOTAL_LATTICE_CELLS]; /* time between hops in each lattice cell */
  int LCells3D[NUM_CELLS_ONESIDE+2][NUM_CELLS_ONESIDE+2][NUM_CELLS_ONESIDE+2]; /* indexes of each lattice cell with padding on top,
                           bottom, left, and right sides to account for the 
                           periodic boundary conditions */
  int Q[TOTAL_LATTICE_CELLS]; /* will hold indexes of the lattice cells organized by their 
               their order in the event queue */
  
  /* filling indexes of lattice cells into an array 
     which is padded on all sides to reflect the
     periodic boundary conditions */
  for(int k=1; k<NUM_CELLS_ONESIDE+1; k++){
    for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
      for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
        LCells3D[i][j][k] = (k-1)*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + (i-1)*NUM_CELLS_ONESIDE + j-1;
      }
    }
  }

  for(int k=1; k<NUM_CELLS_ONESIDE+1; k++){
    for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
      LCells3D[0][j][k] = LCells3D[NUM_CELLS_ONESIDE][j][k];
      LCells3D[NUM_CELLS_ONESIDE+1][j][k] = LCells3D[1][j][k];
    }
  }

  for(int k=1; k<NUM_CELLS_ONESIDE+1; k++){
    for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
      LCells3D[i][0][k] = LCells3D[i][NUM_CELLS_ONESIDE][k];
      LCells3D[i][NUM_CELLS_ONESIDE+1][k] = LCells3D[i][1][k];
    }
  }

  for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
    for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
      LCells3D[i][j][0] = LCells3D[i][j][NUM_CELLS_ONESIDE];
      LCells3D[i][j][NUM_CELLS_ONESIDE+1] = LCells3D[i][j][1];
    }
  }

  /* fill corners with INFINITY since they are irrelevant */
  for(int k=0; k<NUM_CELLS_ONESIDE+2; k++){
    LCells3D[0][0][k] = INFINITY;
    LCells3D[0][NUM_CELLS_ONESIDE+1][k] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][0][k] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][NUM_CELLS_ONESIDE+1][k] = INFINITY;
  }
  for(int j=0; j<NUM_CELLS_ONESIDE+2; j++){
    LCells3D[0][j][0] = INFINITY;
    LCells3D[0][j][NUM_CELLS_ONESIDE+1] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][j][0] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][j][NUM_CELLS_ONESIDE+1] = INFINITY;
  }
  for(int i=0; i<NUM_CELLS_ONESIDE+2; i++){
    LCells3D[i][0][0] = INFINITY;
    LCells3D[i][0][NUM_CELLS_ONESIDE+1] = INFINITY;
    LCells3D[i][NUM_CELLS_ONESIDE+1][0] = INFINITY;
    LCells3D[i][NUM_CELLS_ONESIDE+1][NUM_CELLS_ONESIDE+1] = INFINITY;
  }

  /* Here we fill in the neighbor list of each lattice cell */
  for(int k=0; k<NUM_CELLS_ONESIDE; k++){
    for(int i=0; i<NUM_CELLS_ONESIDE; i++){
      for(int j=0; j<NUM_CELLS_ONESIDE; j++){
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][0] = LCells3D[i][j+1][k+1]; /* North */
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][1] = LCells3D[i+1][j+2][k+1]; /* East */
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][2] = LCells3D[i+2][j+1][k+1]; /* South */
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][3] = LCells3D[i+1][j][k+1]; /* West */
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][4] = LCells3D[i+1][j+1][k+2]; /* Above */
        nn[k*NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE + i*NUM_CELLS_ONESIDE + j][5] = LCells3D[i+1][j+1][k]; /* Below */
      }
    }
  }

// printf("0\t0\t0\t0\t0\n");
/*Loop for running the simulation many times */
for (int iter=0; iter < ITER_MAX; iter++){
  x = 0;
  y = 0;
  z = 0;

    /* Assign the time between hops in each lattice cell.
     Our nanodomain is square with top left corner at
     lattice cell i = 23 + 20*NUM_CELLS_ONESIDE */
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){ /* First we will set all times to TAU_OUT */
    Q[i] = i; /* to be sorted by their order of event occurrence later */
  }

  /* We initialize a uniform distribution of emerin monomer
     proteins over the lattice cells and calculate
     the transition rates for the subvolume KMC method. 
     We then assign a random event time sampled from the 
     probability distribution (master equation). We keep track
     of the total number of proteins Ntot to make sure its constant */ 
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      N[i] =  0;
  }
  N[TOTAL_LATTICE_CELLS/2 + NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE/2 + NUM_CELLS_ONESIDE/2] = TOTAL_BACTERIA_NUM;

  // srand(time(NULL));
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      u1 = (double)rand() / RAND_MAX;
      u2 = (double)rand() / RAND_MAX;
      u3 = (double)rand() / RAND_MAX;
      tt = -TAU_TRAPPED*log(u1);
    //   tt = pow(TAU_TRAPPED*u1,-2.0/3.0);
      th = -TAU_HOP*log(u2);
      W[i] = (double)N[i]/(tt+th)/(2*DIM);
      // W[i] = (double)N[i]/(tt)/4;
      // t[i] = -log(u3)/W[i];
      t[i] = 1.0/W[i];
      Ntot = Ntot + N[i];
  }

  build_min_heap(Q, t, TOTAL_LATTICE_CELLS); /* create priority queue as a binary min heap 
                              organized by the time of the event in each
                              lattice cell */
  int pt = 0; /* event id */
  while(t[0] <= TIME_MAX){ /* we will run the simulation for 4 hours */
    /* assign random integer to picture direction to hop in */
    u4 = rand();
    u5 = (double)rand() / RAND_MAX;
    directionID = u4%(2*DIM);
    hop = (int)ceil(-HOP_AVG*log(u5));
    // if(hop < 1){
    //   hop = 1;
    // }
    // printf("%d\t%d\t%d\t%d\t%f\n",pt+1,x,y,z,t[0]);

    /* grab the index of the lattice cell where the next event will occur */
    lambda = Q[0];
    gamma = lambda;

    switch(directionID){
      case 0 :
        y = y+hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,hop,0,t[0],lambda);
        break;
      case 1 :
        x = x+hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,hop,0,0,t[0],lambda);
        break;
      case 2 :
        y = y-hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,-hop,0,t[0],lambda);
        break;
      case 3 :
        x = x-hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,-hop,0,0,t[0],lambda);
        break;
      case 4 :
        z = z+hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,0,hop,t[0],lambda);
        break;
      case 5 :
        z = z-hop;
        // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,0,-hop,t[0],lambda);
        break;
    }

    /* grab the index of the lattice cell the protein hops into */
    for(int i = 0; i <hop; i++){
      gamma = nn[gamma][directionID];
    }
    printf("%d\t%d\t%f\n",lambda,gamma,t[0]);
    if(pt%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                               where we print a snapshot of the system */
      /* we keep track of the concentration of proteins in the nanodomain
         Nin/Ntot as a function of simulation time t[0] */
      // printf("%d\t%f\t%d\n",pt,t[0],hop);
    }
    /* 1 protein hops out of lattice cell lambda */
    N[lambda] = N[lambda] - 1;
    /* update transition rate in lattice cell lambda */
    u6 = (double)rand() / RAND_MAX;
    u7 = (double)rand() / RAND_MAX;
    tt = -TAU_TRAPPED*log(u6);
    // tt = pow(TAU_TRAPPED*u6,-2.0/3.0);
    th = -TAU_HOP*log(u7);
    W[lambda] = (double)N[lambda]/(tt + th)/(2*DIM);
    // W[lambda] = (double)N[lambda]/(tt)/4;
    /* save the current event time */
    tlambda = t[0];
    /* 1 protein hops into lattice cell gamma */
    N[gamma] = N[gamma] + 1;
    /* update transition rate in lattice cell gamma */
    u8 = (double)rand() / RAND_MAX;
    u9 = (double)rand() / RAND_MAX;
    tt = -TAU_TRAPPED*log(u8);
    // tt = pow(TAU_TRAPPED*u8,-2.0/3.0);
    th = -TAU_HOP*log(u9);
    W[gamma] = (double)N[gamma]/(tt+th)/(2*DIM);
    // W[gamma] = (double)N[gamma]/(tt)/4;

    /* assign random number to calculate new event time in lattice cell lambda */
    u10 = (double)rand() / RAND_MAX;
    // t[0] = tlambda - log(u10) / W[lambda];
    t[0] = tlambda + 1.0/W[lambda];

    /* update order of the priority queue by percolating down the heap*/
    min_heapify(Q, t, 0, TOTAL_LATTICE_CELLS);

    /* search priority queue lattice cell gamma */
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      /* when lattice cell gamma is found in Q
         we update its event time */
      if(Q[i] != gamma){ 
      }
      else{
        /* save event time in lattice cell gamma for comparison */
        tgamma = t[i];
        /* assign random number to calculate new event time in lattice cell gamma */
        u11 = (double)rand() / RAND_MAX;
        // t[i] = tlambda - log(u11)/W[gamma];
        t[i] = tlambda + 1.0/W[gamma];
        /* if new event time is less than old event time
           we use decrease to percolate up the heap */
        if(t[i] < tgamma) {
          decrease_key(Q, t, i);
        }
        else if(t[i] > tgamma) {/* if new event time is greater than old event time
                                   we use heapify to percolate down the heap */
          min_heapify(Q, t, i, TOTAL_LATTICE_CELLS);
        }
        break; /* break out of for loop once lattice cell gamma is found
                  and priority queue is updated */
      }
    }
    Ntot = 0;
    /* recount the number of protein inside the nanodomain and in total */
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      Ntot = Ntot + N[i];
    }
    pt = pt+1; /* increment event id */
  } /* end of while loop after TIME_MAX */
  D = (double)(x*x + y*y + z*z)/(2*DIM)/tlambda/TOTAL_BACTERIA_NUM/(SCALE*SCALE);
  //  printf("\n");
  //  printf("%f\n",D);
  Dsum = Dsum + D;
  Davg = Dsum/(iter+1);
  if(iter%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                               where we print a snapshot of the system */
      /* we keep track of the concentration of proteins in the nanodomain
         Nin/Ntot as a function of simulation time t[0] */
      // printf("%d\t%f\t%d\n",pt,t[0],hop);
    // printf("%d\t%f\n",iter,Davg);
  }
}
// printf("%f",Davg);
printf("%d\t%d\t%f\n",TOTAL_LATTICE_CELLS,SCALE,TIME_MAX);

  
  return 0;
}
