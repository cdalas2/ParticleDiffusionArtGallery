/*******************************************************************************
Subvolume KMC Simulations of freely diffusing monomer emerin proteins in
a domain with a nanodomain of a different diffusion rate due to the clustering
of the emerin proteins there

USAGE

%cc -O3 SKMC.c -o SKMC
%./SKMC > SKMC.out
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "minHeap.h"

#define SNAPSHOT_RATE 1 /* we take a snapshot after every SNAPSHOT_RATE events */
#define TOTAL_PROTEIN_NUM 1600 /* Total number of emerin monomer proteins */
#define NUM_CELLS_ONESIDE 40 /* number of lattice cells along one side of square system domain */
#define TOTAL_LATTICE_CELLS (NUM_CELLS_ONESIDE*NUM_CELLS_ONESIDE) /* Total number of lattice cells */
#define N0 (TOTAL_PROTEIN_NUM/TOTAL_LATTICE_CELLS) /* Initial number of proteins in each lattice */
// #define D_IN 250.0 /* Diffusion coefficient inside nanodomain (nanometers^2/W) */
// #define D_OUT 3500.0 /* Diffusion coefficient outside nanodomain (nanometers^2/W) */
#define D_IN 180 /* Diffusion coefficient inside nanodomain (nanometers^2/W) */
#define D_OUT 540 /* Diffusion coefficient outside nanodomain (nanometers^2/W) */

#define A_IN 12100.0 /* Area of nanodomain (nm^2) */
#define NCPS 11 /* number of lattice cells per side of square nanodomain */
#define A (400.0*400.0) /* Total area of system domain (nm^2) */
#define A_OUT (A - A_IN) /* Area outside nanodomain (nm^2) */
#define LATTICE_CELL_LENGTH 10.0 /* lattice cell length (nm)*/
#define TAU_IN (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_IN) /* time between hops inside nanodomain (W) */
#define TAU_OUT (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_OUT) /* time between hops outside nanodomain (W) */
#define TIME_MAX 600.0
#define SCALE 1 
#define DIM 2
#define LETTERS_CELL_NUM 146

int main() {
  /* directionID will hold lattice neighbor index protein hops into */
  /* lambda will hold the index of the lattice site index it hops out of */
  /* gamma will hold the index of the lattice site index it hops into */
  /* u1,u2,u3,u4 are to hold randoms numbers */
  int directionID, lambda, gamma, u2;
  double u1, u3, u4;

  /* Nin is for the number of proteins inside the nanodomain
     and Ntot is for the total number of proteins across 
     all lattice cells */
  int Nin = 0, Ntot = 0;

  /* tlambda will hold the sampled event time in lattice lambda */
  /* tgamma will hold the sampled event time in lattice gamma */
  double tlambda = 0.0, tgamma = 0.0;

  int N[TOTAL_LATTICE_CELLS]; /* number of emerin monomer proteins in each lattice cell */
  int nn[TOTAL_LATTICE_CELLS][4]; /* nearest neighbors of each lattice cell, we use 
                   periodic boundary conditions at system boundary */
  double W[TOTAL_LATTICE_CELLS]; /* transition rates (diffusion intensities) in each lattice cell */
  double t[TOTAL_LATTICE_CELLS]; /* next queued (sampled) event time in each lattice cell */
  double tau[TOTAL_LATTICE_CELLS]; /* time between hops in each lattice cell */
  int LCells2D[NUM_CELLS_ONESIDE+2][NUM_CELLS_ONESIDE+2]; /* indexes of each lattice cell with padding on top,
                           bottom, left, and right sides to account for the 
                           periodic boundary conditions */
  int Q[TOTAL_LATTICE_CELLS]; /* will hold indexes of the lattice cells organized by their 
               their order in the event queue */
  int Letters[LETTERS_CELL_NUM];

  Letters[0] = 369; Letters[1] =373; Letters[2] =377; Letters[3] =378; Letters[4] =379; Letters[5] =380; Letters[6] =381; 
  Letters[7] =385; Letters[8] = 386; Letters[9] =387; Letters[10] =388; Letters[11] =389;
  Letters[12] =409; Letters[13] =413; Letters[14] =417; Letters[15] =421; Letters[16] =425; Letters[17] =429;
  Letters[18] =449; Letters[19] =453; Letters[20] =457; Letters[21] =465; Letters[22] =469;
  Letters[23] =489; Letters[24] =493; Letters[25] =497; Letters[26] =505;
  Letters[27] =529; Letters[28] =533; Letters[29] =537; Letters[30] =538; Letters[31] =539; Letters[32] =540; Letters[33] =541;
  Letters[34] = 545;
  Letters[35] =569; Letters[36] =573; Letters[37] =581; Letters[38] =585;
  Letters[39] =609; Letters[40] =613; Letters[41] =621; Letters[42] =625; Letters[43] =629;
  Letters[44] =649; Letters[45] =653; Letters[46] =657; Letters[47] =661; Letters[48] =665; Letters[49] =669;
  Letters[50] =689; Letters[51] =690; Letters[52] =691; Letters[53] =692; Letters[54] =693; Letters[55] =697;
  Letters[56] = 698; Letters[57] =699; Letters[58] =700; Letters[59] =701; Letters[60] =705; Letters[61] =706;
  Letters[62] = 707; Letters[63] =708; Letters[64] =709;

  Letters[65] =845; Letters[66] =846; Letters[67] =847; Letters[68] =848; Letters[69] =849; Letters[70] =853;
  Letters[71] = 857; Letters[72] =861; Letters[73] =865; Letters[74] =869; Letters[75] =870; Letters[76] =871;
  Letters[77] =872; Letters[78] =873; 
  Letters[79] =885; Letters[80] =889; Letters[81] =893; Letters[82] =897; Letters[83] =901; Letters[84] =905; Letters[85] =909;
  Letters[144] =913;
  Letters[86] =925; Letters[87] =929; Letters[88] =933; Letters[89] =937; Letters[90] =941; Letters[91] =945; Letters[92] =949;
  Letters[93] =965; Letters[94] =969; Letters[95] =973; Letters[96] =977; Letters[97] =981; Letters[98] =985; Letters[99] =989;
  Letters[100] =1005; Letters[101] =1006; Letters[102] =1007; Letters[103] =1008; Letters[104] =1009; Letters[105] =1013;
  Letters[106] =1014; Letters[107] =1015; Letters[108] =1016; Letters[109] =1017; Letters[110] =1021; Letters[111] =1022;
  Letters[112] = 1023; Letters[113] =1024; Letters[114] =1025; Letters[115] =1029; Letters[116] =1030; Letters[117] =1031; 
  Letters[118] =1032; Letters[119] =1033;
  Letters[120] =1045; Letters[121] =1053; Letters[122] =1057; Letters[123] =1063; Letters[124] =1073;
  Letters[125] =1085; Letters[126] =1093; Letters[127] =1097; Letters[128] =1103; Letters[129] =1113;
  Letters[130] =1125; Letters[131] =1133; Letters[132] =1137; Letters[133] =1143; 
  Letters[145] =1149;
  Letters[134] =1153;
  Letters[135] =1165; Letters[136] =1173; Letters[137] =1177; Letters[138] =1183; Letters[139] =1189;
  Letters[140] = 1190; Letters[141] =1191; Letters[142] =1192; Letters[143] =1193; 
  
  /* Assign the time between hops in each lattice cell.
     Our nanodomain is square with top left corner at
     lattice cell i = 23 + 20*NUM_CELLS_ONESIDE */
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){ /* First we will set all times to TAU_OUT */
    tau[i] = TAU_OUT;
    Q[i] = i; /* to be sorted by their order of event occurrence later */
  }
 
  for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
      for(int j=(23+20*NUM_CELLS_ONESIDE); j<(23+20*NUM_CELLS_ONESIDE+NCPS); j++){
          tau[j + i*NUM_CELLS_ONESIDE] = TAU_IN;      
      }
  }
  //   for(int i=0; i<LETTERS_CELL_NUM; i++){ /* Now we switch in hop times inside nanodomain */
  //     tau[Letters[i]] = TAU_IN;
  // }
  
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
for(int i = 0; i < TOTAL_LATTICE_CELLS; i++){
  N[i] = N0;
}
// N[20] = TOTAL_PROTEIN_NUM;
// for(int i = 0; i < TOTAL_PROTEIN_NUM; i++){
//   N[i] = 1;
// }
  /* We initialize a uniform distribution of emerin monomer
     proteins over the lattice cells and calculate
     the transition rates for the subvolume KMC method. 
     We then assign a random event time sampled from the 
     probability distribution (master equation). We keep track
     of the total number of proteins Ntot to make sure its constant */
  srand(time(NULL));
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      W[i] = (double)N[i]/tau[i]/4; 
      u1 = (double)rand() / RAND_MAX;
      t[i] = -log(u1) / W[i];
      Ntot = Ntot + N[i];
  }

  /* Now we sum up the number of proteins in the nanodomain
     and calculate the total number of proteins outside the
     nanodomain */
  // for(int i=0; i<LETTERS_CELL_NUM; i++) {
  //   Nin = Nin + N[Letters[i]];
  // }

  build_min_heap(Q, t, TOTAL_LATTICE_CELLS); /* create priority queue as a binary min heap 
                              organized by the time of the event in each
                              lattice cell */
  int pt = 0; /* event id */
  while(t[0] <= TIME_MAX){ /* we will run the simulation for 4 hours */
    /* assign random integer to picture direction to hop in */
    u2 = rand();
    directionID = u2%4;

    /* grab the index of the lattice cell where the next event will occur */
    lambda = Q[0];
    /* grab the index of the lattice cell the protein hops into */
    gamma = nn[lambda][directionID];

    if(pt%SNAPSHOT_RATE == 0){ /* this is the snapshot conditional statement
                               where we print a snapshot of the system */
      /* we keep track of the concentration of proteins in the nanodomain
         Nin/Ntot as a function of simulation time t[0] */
      printf("%d\t%d\t%f\t%f\n",lambda,gamma,(double)Nin/Ntot,t[0]);
    }
    /* 1 protein hops out of lattice cell lambda */
    N[lambda] = N[lambda] - 1;
    /* update transition rate in lattice cell lambda */
    W[lambda] = (double)N[lambda]/tau[lambda]/4;
    /* save the current event time */
    tlambda = t[0];
    /* 1 protein hops into lattice cell gamma */
    N[gamma] = N[gamma] + 1;
    /* update transition rate in lattice cell gamma */
    W[gamma] = (double)N[gamma]/tau[gamma]/4;

    /* assign random number to calculate new event time in lattice cell lambda */
    u3 = (double)rand() / RAND_MAX;
    t[0] = tlambda - log(u3) / W[lambda];

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
        u4 = (double)rand() / RAND_MAX;
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

    /* recount the number of protein inside the nanodomain and in total */
    Nin = 0;
    Ntot = 0;
    for(int i=0; i<TOTAL_LATTICE_CELLS; i++){
      Ntot = Ntot + N[i];
    }
    for(int i=0; i<NCPS; i++){
      for(int j=(23+20*NUM_CELLS_ONESIDE); j<(23+20*NUM_CELLS_ONESIDE+NCPS); j++){
          Nin = Nin + N[j + i*NUM_CELLS_ONESIDE];        
      }
    } 
    // for(int i=0; i<LETTERS_CELL_NUM; i++){
    //   Nin = Nin + N[Letters[i]];
    // } 
    
    pt = pt+1; /* increment event id */
  } /* end of while loop after TIME_MAX */
   printf("%d\t%d\t%f\t%f\n",TOTAL_LATTICE_CELLS,SCALE,0.0,TIME_MAX);
  
  return 0;
}
