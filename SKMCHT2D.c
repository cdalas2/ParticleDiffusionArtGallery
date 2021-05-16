/*******************************************************************************
Subvolume KMC Simulations of freely diffusing bacteria in a porous media which
causes it to exhibit a two state motion of hopping and trapping. This code would
be for bacteria in a thin porous film (2D).

USAGE

%cc -O3 SKMCHT2D.c -o SKMCHT2D
%./SKMCHT2D > SKMCHT2D.out
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "minHeap.h"

#define SNAPSHOT_RATE 1000 /* we take a snapshot after every SNAPSHOT_RATE events */
#define TOTAL_BACTERIA_NUM 1 /* Total number of emerin monomer proteins */
#define NUM_CELLS_ONESIDE 5 /* number of lattice cells along one side of square system domain */
#define TOTAL_LATTICE_CELLS (NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE) /* Total number of lattice cells */
#define N0 (TOTAL_BACTERIA_NUM/TOTAL_LATTICE_CELLS) /* Initial number of proteins in each lattice */
#define A (40.0*40.0) /* Total area of system domain (LATTICE_CELL_LENGTH^2) */
#define LATTICE_CELL_LENGTH 1.0 /* lattice cell length (LATTICE_CELL_LENGTH) (=1 micrometer)*/
#define TAU_TRAPPED (0.24) /* time between hops inside nanodomain (s) */
#define TAU_HOP 0.0 /* time between hops outside nanodomain (s) */
#define SCALE 100
#define COORDF 1.0
#define HOP_AVG (3.24*SCALE)
#define TIME_MAX 200.0
#define DIM 2

int main() {
  /* directionID will hold lattice neighbor index protein hops into */
  /* lambda will hold the index of the lattice site index it hops out of */
  /* gamma will hold the index of the lattice site index it hops into */
  /* u1,u2,u3,u4 are to hold randoms numbers */
  int directionID, lambda, gamma, u4, hop, y, x;
  double u1, u2, u3, u5, u6, u7, u8, u9, u10, u11, tt, th, D;

  /* Nin is for the number of proteins inside the nanodomain
     and Ntot is for the total number of proteins across 
     all lattice cells */
  int Ntot = 0;

  /* tlambda will hold the sampled event time in lattice lambda */
  /* tgamma will hold the sampled event time in lattice gamma */
  double tlambda = 0.0, tgamma = 0.0, DAVG = 0.0, DSUM = 0.0;

  int N[TOTAL_LATTICE_CELLS]; /* number of emerin monomer proteins in each lattice cell */
  int nn[TOTAL_LATTICE_CELLS][2*DIM]; /* nearest neighbors of each lattice cell, we use 
                   periodic boundary conditions at system boundary */
  double W[TOTAL_LATTICE_CELLS]; /* transition rates (diffusion intensities) in each lattice cell */
  double t[TOTAL_LATTICE_CELLS]; /* next queued (sampled) event time in each lattice cell */
  double tau[TOTAL_LATTICE_CELLS]; /* time between hops in each lattice cell */
  int LCells2D[NUM_CELLS_ONESIDE+2][NUM_CELLS_ONESIDE+2]; /* indexes of each lattice cell with padding on top,
                           bottom, left, and right sides to account for the 
                           periodic boundary conditions */
  int Q[TOTAL_LATTICE_CELLS]; /* will hold indexes of the lattice cells organized by their 
               their order in the event queue */
  
  /* filling indexes of lattice cells into an array 
     which is padded on all sides to reflect the
     periodic boundary conditions */
  for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
    for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
        LCells2D[i][j] = (i-1)*NUM_CELLS_ONESIDE + j-1;
    }
  }

  for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
      LCells2D[0][j] = j-1 + (NUM_CELLS_ONESIDE-1)*NUM_CELLS_ONESIDE;
      LCells2D[NUM_CELLS_ONESIDE+1][j] = j-1;
  }

  for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
      LCells2D[i][0] = LCells2D[i][NUM_CELLS_ONESIDE];
      LCells2D[i][NUM_CELLS_ONESIDE+1] = LCells2D[i][1];
  }

  /* fill corners with INFINITY since they are irrelevant */
  LCells2D[0][0] = INFINITY;
  LCells2D[0][NUM_CELLS_ONESIDE+1] = INFINITY;
  LCells2D[NUM_CELLS_ONESIDE+1][0] = INFINITY;
  LCells2D[NUM_CELLS_ONESIDE+1][NUM_CELLS_ONESIDE+1] = INFINITY;

  /* Here we fill in the neighbor list of each lattice cell */
  for(int i=0; i<NUM_CELLS_ONESIDE; i++){
    for(int j=0; j<NUM_CELLS_ONESIDE; j++){
        nn[i*NUM_CELLS_ONESIDE + j][0] = LCells2D[i][j+1]; /* up */
        nn[i*NUM_CELLS_ONESIDE + j][1] = LCells2D[i+1][j+2]; /* right */
        nn[i*NUM_CELLS_ONESIDE + j][2] = LCells2D[i+2][j+1]; /* down */
        nn[i*NUM_CELLS_ONESIDE + j][3] = LCells2D[i+1][j]; /* left */
    }
  }

  for (int iter=0; iter < 10000000; iter++){
    x = 0;
    y = 0;

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
    srand(time(NULL));
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
        N[i] =  0;
    }
    N[TOTAL_LATTICE_CELLS/2 + NUM_CELLS_ONESIDE] = TOTAL_BACTERIA_NUM;
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;
        u3 = (double)rand() / RAND_MAX;
        tt = -TAU_TRAPPED*log(u1);
        // tt = pow(1.0*u1,-2.0/3.0);
        th = -TAU_HOP*log(u2);
        W[i] = (double)N[i]/(tt+th)/(2*DIM);
        // W[i] = (double)N[i]/(tt)/4;
        // t[i] = -log(u3)/W[i];
        t[i] = COORDF/W[i];
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
      switch(directionID){
        case 0 :
          y = y+hop;
          break;
        case 1 :
          x = x+hop;
          break;
        case 2 :
          y = y-hop;
          break;
        case 3 :
          x = x-hop;
          break;
      }

      /* grab the index of the lattice cell where the next event will occur */
      lambda = Q[0];
      gamma = lambda;
      /* grab the index of the lattice cell the protein hops into */
      for(int i = 0; i <hop; i++){
        gamma = nn[gamma][directionID];
      }

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
      // tt = pow(1.0*u6,-2.0/3.0);
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
      // tt = pow(1.0*u8,-2.0/3.0);
      th = -TAU_HOP*log(u9);
      W[gamma] = (double)N[gamma]/(tt+th)/4;
      // W[gamma] = (double)N[gamma]/(tt)/4;

      /* assign random number to calculate new event time in lattice cell lambda */
      u10 = (double)rand() / RAND_MAX;
      // t[0] = tlambda - log(u10) / W[lambda];
      t[0] = tlambda + COORDF/W[lambda];

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
          t[i] = tlambda + COORDF/W[gamma];
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
    D = (double)(x*x + y*y)/(2*DIM)/tlambda/TOTAL_BACTERIA_NUM/(SCALE*SCALE);
    //  printf("\n");
    //  printf("%f\n",D);
    DSUM = DSUM + D;
    DAVG = DSUM/(iter+1);
    if(iter%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                                where we print a snapshot of the system */
        /* we keep track of the concentration of proteins in the nanodomain
          Nin/Ntot as a function of simulation time t[0] */
        // printf("%d\t%f\t%d\n",pt,t[0],hop);
      printf("%d\t%f\n",iter,DAVG);
    }
  }  

  return 0;
}
