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
#define D_IN 30.0 /* Diffusion coefficient inside nanodomain (nanometers^2/W) */
#define D_OUT 1000.0 /* Diffusion coefficient outside nanodomain (nanometers^2/W) */

#define A_IN 12100.0 /* Area of nanodomain (nm^2) */
#define NCPS 11 /* number of lattice cells per side of square nanodomain */
#define A (400.0*400.0) /* Total area of system domain (nm^2) */
#define A_OUT (A - A_IN) /* Area outside nanodomain (nm^2) */
#define LATTICE_CELL_LENGTH 10.0 /* lattice cell length (nm)*/
#define TAU_IN (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_IN) /* time between hops inside nanodomain (W) */
#define TAU_OUT (LATTICE_CELL_LENGTH*LATTICE_CELL_LENGTH/D_OUT) /* time between hops outside nanodomain (W) */
#define TIME_MAX 14500.0
#define STEP_MAX 1000000.0
#define SCALE 1 
#define DIM 2
#define HEAD_OUTLINE_CELL_NUM 272
#define HEAD_INSIDE_CELL_NUM 162
#define CHAIN_CELL_NUM 136
#define BILAYER_EMPTY_SPACE_INSIDE 292
#define BILAYER_ROW_SPACING 16

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
  int SubDomainCellNum = 0*HEAD_OUTLINE_CELL_NUM + 1*HEAD_INSIDE_CELL_NUM + 0*CHAIN_CELL_NUM + 0*BILAYER_EMPTY_SPACE_INSIDE;
  int SubDomainCells[SubDomainCellNum];

  int offset = 40*9;

// //TOP HEAD OUTLINE
//   SubDomainCells[0] = 0; SubDomainCells[1] =3; SubDomainCells[2] =4; SubDomainCells[3] =5; SubDomainCells[4] =6; SubDomainCells[5] =9; SubDomainCells[6] =10; 
//   SubDomainCells[7] =11; SubDomainCells[8] = 12; SubDomainCells[9] =15; SubDomainCells[10] =16; SubDomainCells[11] =17;
//   SubDomainCells[12] =18; SubDomainCells[13] =21; SubDomainCells[14] =22; SubDomainCells[15] =23; SubDomainCells[16] =24; SubDomainCells[17] =27;
//   SubDomainCells[18] =28; SubDomainCells[19] =29; SubDomainCells[20] =30; SubDomainCells[21] =33; SubDomainCells[22] =34;
//   SubDomainCells[23] =35; SubDomainCells[24] =36; SubDomainCells[25] =39; SubDomainCells[26] =40;
//   SubDomainCells[27] =41; SubDomainCells[28] =42; SubDomainCells[29] =43; SubDomainCells[30] =46; SubDomainCells[31] =47; SubDomainCells[32] =48; SubDomainCells[33] =49;
//   SubDomainCells[34] = 52;
//   SubDomainCells[35] =53; SubDomainCells[36] =54; SubDomainCells[37] =55; SubDomainCells[38] =58;
//   SubDomainCells[39] =59; SubDomainCells[40] =60; SubDomainCells[41] =61; SubDomainCells[42] =64; SubDomainCells[43] =65;
//   SubDomainCells[44] =66; SubDomainCells[45] =67; SubDomainCells[46] =70; SubDomainCells[47] =71; SubDomainCells[48] =72; SubDomainCells[49] =73;
//   SubDomainCells[50] =76; SubDomainCells[51] =77; SubDomainCells[52] =78; SubDomainCells[53] =79; SubDomainCells[54] =81; SubDomainCells[55] =82;
//   SubDomainCells[56] = 87; SubDomainCells[57] =88; SubDomainCells[58] =93; SubDomainCells[59] =94; SubDomainCells[60] =99; SubDomainCells[61] =100;
//   SubDomainCells[62] = 105; SubDomainCells[63] =106; SubDomainCells[64] =111;

//   SubDomainCells[65] =112; SubDomainCells[66] =117; SubDomainCells[67] =118; SubDomainCells[68] =121; SubDomainCells[69] =122; SubDomainCells[70] =127;
//   SubDomainCells[71] = 128; SubDomainCells[72] =133; SubDomainCells[73] =134; SubDomainCells[74] =139; SubDomainCells[75] =140; SubDomainCells[76] =145;
//   SubDomainCells[77] =146; SubDomainCells[78] =151; 
//   SubDomainCells[79] =152; SubDomainCells[80] =157; SubDomainCells[81] =158; SubDomainCells[82] =160; SubDomainCells[83] =161; SubDomainCells[84] =162; SubDomainCells[85] =163;
//   SubDomainCells[86] =166;
//   SubDomainCells[87] =167; SubDomainCells[88] =168; SubDomainCells[89] =169; SubDomainCells[90] =170; SubDomainCells[91] =173; SubDomainCells[92] =174; SubDomainCells[93] =175;
//   SubDomainCells[94] =176; SubDomainCells[95] =179; SubDomainCells[96] =180; SubDomainCells[97] =181; SubDomainCells[98] =182; SubDomainCells[99] =185; SubDomainCells[100] =186;
//   SubDomainCells[101] =187; SubDomainCells[102] =188; SubDomainCells[103] =191; SubDomainCells[104] =192; SubDomainCells[105] =193; SubDomainCells[106] =194;
//   SubDomainCells[107] =197; SubDomainCells[108] =198; SubDomainCells[109] =199; SubDomainCells[110] =200; SubDomainCells[111] =203; SubDomainCells[112] =204;
//   SubDomainCells[113] = 205; SubDomainCells[114] =206; SubDomainCells[115] =209; SubDomainCells[116] =210; SubDomainCells[117] =211; SubDomainCells[118] =212; 
//   SubDomainCells[119] =215; SubDomainCells[120] =216;
//   SubDomainCells[121] =217; SubDomainCells[122] =218; SubDomainCells[123] =221; SubDomainCells[124] =222; SubDomainCells[125] =223;
//   SubDomainCells[126] =224; SubDomainCells[127] =227; SubDomainCells[128] =228; SubDomainCells[129] =229; SubDomainCells[130] =230;
//   SubDomainCells[131] =233; SubDomainCells[132] =234; SubDomainCells[133] =235; SubDomainCells[134] =236; 
//   SubDomainCells[135] =239;

// //BOTTOM HEAD OUTLINE
//   for(int i = 0; i < (HEAD_OUTLINE_CELL_NUM/2); i++){
//       SubDomainCells[(HEAD_OUTLINE_CELL_NUM/2)+i] = SubDomainCells[i] + NUM_CELLS_ONESIDE*BILAYER_ROW_SPACING;
//   }

  //TOP HEAD INSIDE
  SubDomainCells[0] = 44; SubDomainCells[1] =45; SubDomainCells[2] =50; SubDomainCells[3] =51; SubDomainCells[4] =56; SubDomainCells[5] =57; SubDomainCells[6] =62; 
  SubDomainCells[7] =11; SubDomainCells[8] = 63; SubDomainCells[9] =68; SubDomainCells[10] =69; SubDomainCells[11] =74;
  SubDomainCells[12] =18; SubDomainCells[13] =75; SubDomainCells[14] =80; SubDomainCells[15] =83; SubDomainCells[16] =84; SubDomainCells[17] =85;
  SubDomainCells[18] =28; SubDomainCells[19] =86; SubDomainCells[20] =89; SubDomainCells[21] =90; SubDomainCells[22] =91;
  SubDomainCells[23] =35; SubDomainCells[24] =92; SubDomainCells[25] =95; SubDomainCells[26] =96;
  SubDomainCells[27] =41; SubDomainCells[28] =97; SubDomainCells[29] =98; SubDomainCells[30] =101; SubDomainCells[31] =102; SubDomainCells[32] =103; SubDomainCells[33] =104;
  SubDomainCells[34] = 107;
  SubDomainCells[35] =108; SubDomainCells[36] =109; SubDomainCells[37] =110; SubDomainCells[38] =113;
  SubDomainCells[39] =114; SubDomainCells[40] =115; SubDomainCells[41] =116; SubDomainCells[42] =119; SubDomainCells[43] =120;
  SubDomainCells[44] =123; SubDomainCells[45] =124; SubDomainCells[46] =125; SubDomainCells[47] =126; SubDomainCells[48] =129; SubDomainCells[49] =130;
  SubDomainCells[50] =131; SubDomainCells[51] =132; SubDomainCells[52] =135; SubDomainCells[53] =136; SubDomainCells[54] =137; SubDomainCells[55] =138;
  SubDomainCells[56] = 141; SubDomainCells[57] =142; SubDomainCells[58] =143; SubDomainCells[59] =144; SubDomainCells[60] =147; SubDomainCells[61] =148;
  SubDomainCells[62] = 149; SubDomainCells[63] =150; SubDomainCells[64] =153;

  SubDomainCells[65] =154; SubDomainCells[66] =155; SubDomainCells[67] =156; SubDomainCells[68] =159; SubDomainCells[69] =164; SubDomainCells[70] =165;
  SubDomainCells[71] = 170; SubDomainCells[72] =171; SubDomainCells[73] =176; SubDomainCells[74] =177; SubDomainCells[75] =182; SubDomainCells[76] =183;
  SubDomainCells[77] =188; SubDomainCells[78] =189; 
  SubDomainCells[79] =194; SubDomainCells[80] =195;

//BOTTOM HEAD INSIDE
  for(int i = 0; i < (HEAD_INSIDE_CELL_NUM/2); i++){
      SubDomainCells[(HEAD_INSIDE_CELL_NUM/2)+i] = SubDomainCells[i] + NUM_CELLS_ONESIDE*BILAYER_ROW_SPACING;
  }

// //TOP CHAINS
//   SubDomainCells[0] = 240; SubDomainCells[1] =243; SubDomainCells[2] =244; SubDomainCells[3] =245; SubDomainCells[4] =246; SubDomainCells[5] =249; SubDomainCells[6] =250; 
//   SubDomainCells[7] =251; SubDomainCells[8] = 252; SubDomainCells[9] =255; SubDomainCells[10] =256; SubDomainCells[11] =257;
//   SubDomainCells[12] =258; SubDomainCells[13] =261; SubDomainCells[14] =262; SubDomainCells[15] =263; SubDomainCells[16] =264; SubDomainCells[17] =267;
//   SubDomainCells[18] =268; SubDomainCells[19] =269; SubDomainCells[20] =270; SubDomainCells[21] =273; SubDomainCells[22] =274;
//   SubDomainCells[23] =275; SubDomainCells[24] =276; SubDomainCells[25] =279; SubDomainCells[26] =280;
//   SubDomainCells[27] =283; SubDomainCells[28] =286; SubDomainCells[29] =289; SubDomainCells[30] =292; SubDomainCells[31] =295; SubDomainCells[32] =298; SubDomainCells[33] =301;
//   SubDomainCells[34] = 304;
//   SubDomainCells[35] =307; SubDomainCells[36] =310; SubDomainCells[37] =313; SubDomainCells[38] =316;
//   SubDomainCells[39] =319; SubDomainCells[40] =320; SubDomainCells[41] =323; SubDomainCells[42] =326; SubDomainCells[43] =329;
//   SubDomainCells[44] =332; SubDomainCells[45] =335; SubDomainCells[46] =338; SubDomainCells[47] =341; SubDomainCells[48] =344; SubDomainCells[49] =347;
//   SubDomainCells[50] =350; SubDomainCells[51] =353; SubDomainCells[52] =356; SubDomainCells[53] =359; SubDomainCells[54] =360; SubDomainCells[55] =363;
//   SubDomainCells[56] = 366; SubDomainCells[57] =369; SubDomainCells[58] =372; SubDomainCells[59] =375; SubDomainCells[60] =378; SubDomainCells[61] =381;
//   SubDomainCells[62] = 384; SubDomainCells[63] =387; SubDomainCells[64] =390;

//   SubDomainCells[65] =393; SubDomainCells[66] =396; SubDomainCells[67] =399;

// //BOTTOM CHAINS
//   for(int i = 0; i < 14; i++){
//       SubDomainCells[(CHAIN_CELL_NUM/2) + i] = SubDomainCells[i+26+14+14] + 3*NUM_CELLS_ONESIDE;
//       SubDomainCells[(CHAIN_CELL_NUM/2) + i + 14] = SubDomainCells[i+26+14] + 5*NUM_CELLS_ONESIDE;
//       SubDomainCells[(CHAIN_CELL_NUM/2) + i + 14 + 14] = SubDomainCells[i+26] + 7*NUM_CELLS_ONESIDE;
//   }
//   for(int i = 0; i < 26; i++){
//       SubDomainCells[(CHAIN_CELL_NUM-26)+i] = SubDomainCells[i] + 9*NUM_CELLS_ONESIDE;
//   }

//   //BILAYER EMPTY SPACE INSIDE
//   SubDomainCells[0] = 201; SubDomainCells[1] =202; SubDomainCells[2] =207; SubDomainCells[3] =208; SubDomainCells[4] =213; SubDomainCells[5] =214; SubDomainCells[6] =219; 
//   SubDomainCells[7] =220; SubDomainCells[8] = 225; SubDomainCells[9] =226; SubDomainCells[10] =231; SubDomainCells[11] =232;
//   SubDomainCells[12] =237; SubDomainCells[13] =238; SubDomainCells[14] =241; SubDomainCells[15] =242; SubDomainCells[16] =247; SubDomainCells[17] =248;
//   SubDomainCells[18] =253; SubDomainCells[19] =254; SubDomainCells[20] =259; SubDomainCells[21] =260; SubDomainCells[22] =265;
//   SubDomainCells[23] =266; SubDomainCells[24] =271; SubDomainCells[25] =272; SubDomainCells[26] =277;
//   SubDomainCells[27] =278; SubDomainCells[28] =281; SubDomainCells[29] =282; SubDomainCells[30] =284; SubDomainCells[31] =285; SubDomainCells[32] =287; SubDomainCells[33] =288;
//   SubDomainCells[34] = 290;
//   SubDomainCells[35] =291; SubDomainCells[36] =293; SubDomainCells[37] =294; SubDomainCells[38] =296;
//   SubDomainCells[39] =297; SubDomainCells[40] =299; SubDomainCells[41] =300; SubDomainCells[42] =302; SubDomainCells[43] =303;
//   SubDomainCells[44] =305; SubDomainCells[45] =306; SubDomainCells[46] =308; SubDomainCells[47] =309; SubDomainCells[48] =311; SubDomainCells[49] =312;
//   SubDomainCells[50] =314; SubDomainCells[51] =315; SubDomainCells[52] =317; SubDomainCells[53] =318; SubDomainCells[54] =321; SubDomainCells[55] =322;
//   SubDomainCells[56] = 324; SubDomainCells[57] =325; SubDomainCells[58] =327; SubDomainCells[59] =328; SubDomainCells[60] =330; SubDomainCells[61] =331;
//   SubDomainCells[62] = 333; SubDomainCells[63] =334; SubDomainCells[64] =336;

//   SubDomainCells[65] =337; SubDomainCells[66] =339; SubDomainCells[67] =340; SubDomainCells[68] =342; SubDomainCells[69] =343; SubDomainCells[70] =345;
//   SubDomainCells[71] = 346; SubDomainCells[72] =348; SubDomainCells[73] =349; SubDomainCells[74] =351; SubDomainCells[75] =352; SubDomainCells[76] =354;
//   SubDomainCells[77] =355; SubDomainCells[78] =357; 
//   SubDomainCells[79] =358; SubDomainCells[80] =361; SubDomainCells[81] =362; SubDomainCells[82] =364; SubDomainCells[83] =365; SubDomainCells[84] =367; SubDomainCells[85] =368;
//   SubDomainCells[86] =370;
//   SubDomainCells[87] =371; SubDomainCells[88] =373; SubDomainCells[89] =374; SubDomainCells[90] =376; SubDomainCells[91] =377; SubDomainCells[92] =379; SubDomainCells[93] =380;
//   SubDomainCells[94] =382; SubDomainCells[95] =383; SubDomainCells[96] =385; SubDomainCells[97] =386; SubDomainCells[98] =388; SubDomainCells[99] =389; SubDomainCells[100] =391;
//   SubDomainCells[101] =392; SubDomainCells[102] =394; SubDomainCells[103] =395; SubDomainCells[104] =397; SubDomainCells[105] =398; 
//   for(int i = 0; i < 80; i++){
//       SubDomainCells[106+i] = 400 + i;
//   }
//   for(int i = 0; i < 26; i++){
//       SubDomainCells[(BILAYER_EMPTY_SPACE_INSIDE/2) + i + 40] = SubDomainCells[i+14+14+26+26] + 3*NUM_CELLS_ONESIDE;
//       SubDomainCells[(BILAYER_EMPTY_SPACE_INSIDE/2) + i + 40 + 26] = SubDomainCells[i+14+14+26] + 5*NUM_CELLS_ONESIDE;
//       SubDomainCells[(BILAYER_EMPTY_SPACE_INSIDE/2) + i + 40 + 26 + 26] = SubDomainCells[i+14+14] + 7*NUM_CELLS_ONESIDE;
//   }
//   for(int i = 0; i < 14; i++){
//       SubDomainCells[(BILAYER_EMPTY_SPACE_INSIDE/2) + i + 40 + 26 + 26 + 26] = SubDomainCells[i+14] + 9*NUM_CELLS_ONESIDE;
//       SubDomainCells[(BILAYER_EMPTY_SPACE_INSIDE/2) + i + 40 + 26 + 26 + 14] = SubDomainCells[i] + 11*NUM_CELLS_ONESIDE;

//   }


  for(int i =0; i < SubDomainCellNum; i++){
      SubDomainCells[i] = SubDomainCells[i] + offset;
  }
  /* Assign the time between hops in each lattice cell.
     Our nanodomain is square with top left corner at
     lattice cell i = 23 + 20*NUM_CELLS_ONESIDE */
  for(int i=0; i<TOTAL_LATTICE_CELLS; i++){ /* First we will set all times to TAU_OUT */
    tau[i] = TAU_OUT;
    Q[i] = i; /* to be sorted by their order of event occurrence later */
  }
 
//   for(int i=0; i<NCPS; i++){ /* Now we switch in hop times inside nanodomain */
//       for(int j=(23+20*NUM_CELLS_ONESIDE); j<(23+20*NUM_CELLS_ONESIDE+NCPS); j++){
//           tau[j + i*NUM_CELLS_ONESIDE] = TAU_IN;      
//       }
//   }
    for(int i=0; i<SubDomainCellNum; i++){ /* Now we switch in hop times inside nanodomain */
      tau[SubDomainCells[i]] = TAU_IN;
  }
  
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
  N[i] = 0;
}
N[800] = TOTAL_PROTEIN_NUM;
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
  for(int i=0; i<SubDomainCellNum; i++) {
    Nin = Nin + N[SubDomainCells[i]];
  }

  build_min_heap(Q, t, TOTAL_LATTICE_CELLS); /* create priority queue as a binary min heap 
                              organized by the time of the event in each
                              lattice cell */
  int pt = 0; /* event id */
  while(pt <= STEP_MAX){ /* we will run the simulation for 4 hours */
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
    // for(int i=0; i<NCPS; i++){
    //   for(int j=(23+20*NUM_CELLS_ONESIDE); j<(23+20*NUM_CELLS_ONESIDE+NCPS); j++){
    //       Nin = Nin + N[j + i*NUM_CELLS_ONESIDE];        
    //   }
    // } 
    for(int i=0; i<SubDomainCellNum; i++){
      Nin = Nin + N[SubDomainCells[i]];
    } 
    
    pt = pt+1; /* increment event id */
  } /* end of while loop after TIME_MAX */
   printf("%d\t%d\t%f\t%f\n",TOTAL_LATTICE_CELLS,SCALE,0.0,STEP_MAX);
  
  return 0;
}
