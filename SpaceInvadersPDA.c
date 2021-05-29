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
#define NUM_CELLS_ONESIDE 40 /* number of lattice cells along one side of square system domain */
#define TOTAL_LATTICE_CELLS (NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE) /* Total number of lattice cells */
#define TOTAL_BACTERIA_NUM 1600 /* Total number of emerin monomer proteins */
#define N0 (TOTAL_BACTERIA_NUM/TOTAL_LATTICE_CELLS) /* Initial number of proteins in each lattice */
#define NCPS 18 /* number of lattice cells per side of square nanodomain */
// #define V (6.0*6.0*6.0) /* Total area of system domain (LATTICE_CELL_LENGTH^2) */
#define D_IN 180.0
#define D_OUT 540.0
#define LATTICE_CELL_LENGTH 10.0 /* lattice cell length (=10.0 micrometer/SCALE)*/
#define TAU_IN (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_IN)/* time between hops inside nanodomain (s) */
#define TAU_OUT (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_OUT) /* time between hops outside nanodomain (s). we set to zero since none were reported in the experiments. */
#define SCALE 1 /*re-scales the size of the lattice cells */
#define TIME_MAX 600.0 /*simulation time limit */
#define ITER_MAX 1 /*number of iterations of the simulation */

#define SDSPEED 4
#define DIM 2 /* Our domain is 2D */

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
  int directionID, lambda, gamma, u2, y, x, Nr, lr, ls, pt, interval, SDBase, jr;
  double u1, u3, u4, u5, u6, D, tr,tF;

  /* Ntot is for the total number of proteins across 
     all lattice cells */
  int Ntot = 0;
  int Nin = 0;

  /* tlambda will hold the sampled event time in lattice lambda */
  /* tgamma will hold the sampled event time in lattice gamma */
  /* Dsum will hold the diffusivities summed over iterations for averaging in later step */
  /* Davg will hold the diffusivity averaged over iterations */
  double tlambda = 0.0, tgamma = 0.0, Dsum = 0.0, Davg = 0.0;

  int N[TOTAL_LATTICE_CELLS]; /* number of emerin monomer proteins in each lattice cell */
  int nn[TOTAL_LATTICE_CELLS][2*DIM]; /* nearest neighbors of each lattice cell, we use 
                   periodic boundary conditions at system boundary */
  double W[TOTAL_LATTICE_CELLS]; /* transition rates (diffusion intensities) in each lattice cell */
  double t[TOTAL_LATTICE_CELLS]; /* next queued (sampled) event time in each lattice cell */
  double tau[TOTAL_LATTICE_CELLS]; /* time between hops in each lattice cell */
  int LCells3D[NUM_CELLS_ONESIDE+2][NUM_CELLS_ONESIDE+2]; /* indexes of each lattice cell with padding on top,
                           bottom, left, and right sides to account for the 
                           periodic boundary conditions */
  int Q[TOTAL_LATTICE_CELLS]; /* will hold indexes of the lattice cells organized by their 
               their order in the event queue */
  
  /* filling indexes of lattice cells into an array 
     which is padded on all sides to reflect the
     periodic boundary conditions */
    for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
      for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
        LCells3D[i][j] = (i-1)*NUM_CELLS_ONESIDE + j-1;
      }
    }

    for(int j=1; j<NUM_CELLS_ONESIDE+1; j++){
      LCells3D[0][j] = LCells3D[NUM_CELLS_ONESIDE][j];
      LCells3D[NUM_CELLS_ONESIDE+1][j] = LCells3D[1][j];
    }

    for(int i=1; i<NUM_CELLS_ONESIDE+1; i++){
      LCells3D[i][0] = LCells3D[i][NUM_CELLS_ONESIDE];
      LCells3D[i][NUM_CELLS_ONESIDE+1] = LCells3D[i][1];
    }

  /* fill corners with INFINITY since they are irrelevant */
    LCells3D[0][0] = INFINITY;
    LCells3D[0][NUM_CELLS_ONESIDE+1] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][0] = INFINITY;
    LCells3D[NUM_CELLS_ONESIDE+1][NUM_CELLS_ONESIDE+1] = INFINITY;

  /* Here we fill in the neighbor list of each lattice cell */
    for(int i=0; i<NUM_CELLS_ONESIDE; i++){
      for(int j=0; j<NUM_CELLS_ONESIDE; j++){
        nn[i*NUM_CELLS_ONESIDE + j][0] = LCells3D[i][j+1]; /* Up */
        nn[i*NUM_CELLS_ONESIDE + j][1] = LCells3D[i+1][j+2]; /* Right */
        nn[i*NUM_CELLS_ONESIDE + j][2] = LCells3D[i+2][j+1]; /* Down */
        nn[i*NUM_CELLS_ONESIDE + j][3] = LCells3D[i+1][j]; /* Left */
      }
    }

  // printf("0\t0\t0\t0\t0\n");
  /*Loop for running the simulation many times */
  for (int iter=0; iter < ITER_MAX; iter++){
    pt = 0; /* event id */
    interval = 1;
    SDBase = 0;
    tr = 0.0;
    x = 0;
    y = 0;

    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){ /* First we will set all times to TAU_OUT */
      tau[i] = TAU_OUT;
    }
    for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
        for(int j=0; j<NCPS; j++){
            tau[j + i*NUM_CELLS_ONESIDE] = TAU_IN;  
        }
    } 
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
        N[i] =  N0;
    }
    
    // N[0] = TOTAL_BACTERIA_NUM;

    srand(time(NULL));
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
        u1 = (double)rand() / RAND_MAX;
        W[i] = (double)N[i]/tau[i]/(2*DIM);
        t[i] = -log(u1)/W[i];
        Ntot = Ntot + N[i];
    }
    for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
        for(int j=0; j<NCPS; j++){
            Nin = Nin + N[j + i*NUM_CELLS_ONESIDE];  
        }
    } 

    build_min_heap(Q, t, TOTAL_LATTICE_CELLS); /* create priority queue as a binary min heap 
                                organized by the time of the event in each
                                lattice cell */
    tF = 0.0;
    while(t[0] <= TIME_MAX){ /* we will run the simulation for 4 hours */
      /* assign random integer to picture direction to hop in */
      // while(t[0] < tF){
      //   min_heapify(Q, t, 0, TOTAL_LATTICE_CELLS);
      // }

      if(t[0] < (double)interval/SDSPEED){
        u2 = rand();
        directionID = u2%(2*DIM);

        /* grab the index of the lattice cell where the next event will occur */
        lambda = Q[0];
        gamma = nn[lambda][directionID];

        switch(directionID){
          case 0 :
            y = y+1;
            // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,hop,0,t[0],lambda);
            break;
          case 1 :
            x = x+1;
            // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,hop,0,0,t[0],lambda);
            break;
          case 2 :
            y = y-1;
            // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,0,-hop,0,t[0],lambda);
            break;
          case 3 :
            x = x-1;
            // printf("%d\t%d\t%d\t%d\t%f\t%d\n",pt+1,-hop,0,0,t[0],lambda);
            break;
        }
        /* 1 protein hops out of lattice cell lambda */
        N[lambda] = N[lambda] - 1;
        /* update transition rate in lattice cell lambda */
        W[lambda] = (double)N[lambda]/tau[lambda]/(2*DIM);
        /* save the current event time */
        tlambda = t[0];
        /* 1 protein hops into lattice cell gamma */
        N[gamma] = N[gamma] + 1;
        /* update transition rate in lattice cell gamma */
        W[gamma] = (double)N[gamma]/tau[gamma]/(2*DIM);

        /* assign random number to calculate new event time in lattice cell lambda */
        u3 = (double)rand() / RAND_MAX;
        // t[0] = tlambda - log(u10) / W[lambda];
        t[0] = tlambda - log(u3)/W[lambda];

        /* update order of the priority queue by percolating down the heap*/
        min_heapify(Q, t, 0, TOTAL_LATTICE_CELLS);

        /* search for priority queue lattice cell gamma */
        for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
          /* when lattice cell gamma is found in Q
            we update its event time */
          if(Q[i] != gamma){ 
          }
          else{
            /* save event time in lattice cell gamma for comparison */
            tgamma = t[i];
            /* assign random number to calculate new event time in lattice cell gamma */
            u4 = (double)rand() / RAND_MAX;
            // t[i] = tlambda - log(u11)/W[gamma];
            t[i] = tlambda - log(u4)/W[gamma];
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

        if(pt%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                                  where we print a snapshot of the system */
          /* we keep track of the concentration of proteins in the nanodomain
            Nin/Ntot as a function of simulation time t[0] */
          // printf("%d\t%f\t%d\n",pt,t[0],hop);
        
          Ntot = 0;
          Nin = 0;
          for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
            Ntot = Ntot + N[i];
          }
          for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
            jr = SDBase + i*NUM_CELLS_ONESIDE;
            for(int j=0; j<NCPS; j++){
              Nin = Nin + N[jr];
              jr = nn[jr][1];
            }
          } 

          printf("%d\t%d\t%f\n",lambda,gamma,tlambda);
        }
        tF = tlambda;
      }
      else{
          for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
            tau[SDBase + i*NUM_CELLS_ONESIDE] = TAU_OUT; 
            lr = SDBase + i*NUM_CELLS_ONESIDE;
            for(int j=0; j<NCPS; j++){
              lr = nn[lr][1];
            }

            tau[lr] = TAU_IN;  
            ls = nn[lr][3];
            N[lr]  = N[lr] + N[ls];
            W[lr] = (double)N[lr]/tau[lr]/(2*DIM);
            /* assign random number to calculate new event time in lattice cell gamma */
            u5 = (double)rand() / RAND_MAX;
            // t[i] = tlambda - log(u11)/W[gamma];
            tr = (double)interval/SDSPEED - log(u5)/W[lr];
            for(int j=0; j<TOTAL_LATTICE_CELLS; j++){
              if( Q[j] != lr ){ 
              }
              else{
                /* save event time in lattice cell gamma for comparison */
                if(tr < t[j]) {
                t[j] = tr;
                decrease_key(Q, t, j);
                }
                else if(tr > t[j]) {/* if new event time is greater than old event time
                                        we use heapify to percolate down the heap */
                t[j] = tr;
                min_heapify(Q, t, j, TOTAL_LATTICE_CELLS);
                }
                break; /* break out of for loop once lattice cell gamma is found
                        and priority queue is updated */
              }
            }
          }
          for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
            lr = SDBase + i*NUM_CELLS_ONESIDE;
            for(int j=0; j<NCPS; j++){
              lr = nn[lr][1];
            }
            for(int j=1; j<NCPS; j++){         
              lr = nn[lr][3];
              ls = nn[lr][3];

              N[lr]  = N[ls];
              W[lr] = (double)N[lr]/tau[lr]/(2*DIM);
              u6 = (double)rand() / RAND_MAX;
              tr = (double)interval/SDSPEED - log(u6)/W[lr];
              for(int q=0; q<TOTAL_LATTICE_CELLS; q++){
                if( Q[q] != lr ){ 
                }
                else{
                  /* save event time in lattice cell gamma for comparison */
                  if(tr < t[q]) {
                    t[q] = tr;
                    decrease_key(Q, t, q);
                  }
                  else if(tr > t[q]) {/* if new event time is greater than old event time
                                            we use heapify to percolate down the heap */
                    t[q] = tr;
                    min_heapify(Q, t, q, TOTAL_LATTICE_CELLS);
                  }
                  break; /* break out of for loop once lattice cell gamma is found
                            and priority queue is updated */
                }
              }
            }
            lr = SDBase + i*NUM_CELLS_ONESIDE;
            N[lr] = 0;
            W[lr] = 0;
            tr = INFINITY;
            for(int q=0; q<TOTAL_LATTICE_CELLS; q++){
              if( Q[q] != lr ){ 
              }
              else{
                if(tr > t[q]) {
                  t[q] = tr;
                }
                min_heapify(Q, t, q,TOTAL_LATTICE_CELLS);
                break; /* break out of for loop once lattice cell gamma is found
                          and priority queue is updated */
              }
            }
          }
        SDBase = nn[SDBase][1];
        if(pt%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                                  where we print a snapshot of the system */
          /* we keep track of the concentration of proteins in the nanodomain
            Nin/Ntot as a function of simulation time t[0] */
          // printf("%d\t%f\t%d\n",pt,t[0],hop);
        
          Ntot = 0;
          Nin = 0;
          for(int i=0; i<TOTAL_LATTICE_CELLS;i++){
            Ntot = Ntot + N[i];
          }
            for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
              jr = SDBase + i*NUM_CELLS_ONESIDE;
              for(int j=0; j<NCPS; j++){
                Nin = Nin + N[jr];
                jr = nn[jr][1];
              }
            } 
          printf("%d\t%d\t%f\n",10*TOTAL_LATTICE_CELLS,10*TOTAL_LATTICE_CELLS,(double)interval/SDSPEED);
        }
        tF = (double)interval/SDSPEED;
        interval = interval + 1;
      }
      pt = pt+1; /* increment event id */
    } /* end of while loop after TIME_MAX */
    D = (double)(x*x + y*y)/(2*DIM)/tlambda/TOTAL_BACTERIA_NUM/(SCALE*SCALE);
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
  printf("%d\t%d\t%f\n",TOTAL_LATTICE_CELLS,SDSPEED,TIME_MAX);

  
  return 0;
}
