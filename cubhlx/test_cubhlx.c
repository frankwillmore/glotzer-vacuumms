#include <stdio.h>
#include <math.h>

//      SUBROUTINE CUBHLX(START,ROTS,HUE,GAMMA,
      
//      number of levels

int main(){
/*
  INTEGER   NLEV,I,NLO,NHI
  REAL      START,ROTS,HUE,GAMMA
  REAL      RED(NLEV),GRN(NLEV),BLU(NLEV)
  REAL      PI,FRACT,ANGLE,AMP
*/
  int       NLO,NHI, NLEV=256;
  float     START=3,ROTS=1,HUE=1,GAMMA=0.333;
  float     RED[NLEV],GRN[NLEV],BLU[NLEV];

  cubhlx_(&START,&ROTS,&HUE,&GAMMA,&NLEV,&RED,&GRN,&BLU,&NLO,&NHI);

  int i;

  for (i=0; i<NLEV; i++){
    printf("red %f\n", RED[i]);
  }

  return 0;
}
