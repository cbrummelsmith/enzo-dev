/***********************************************************************
/
/  GRID: ADD Stellar Feedback to Elliptical Galaxies (Hernquist profile)
/
/  written by: Yuan Li
/  date:       Aug, 2015
/  modified1: 
/
/  PURPOSE: Stellar wind adds mass, and SN Ia adds heat
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::OldStarFeedback()
{
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
                                             Vel3Num, TENum) == FAIL)   ///this or thisgrid
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.");


//ClusterSMBHBCG is M_* here
//EllipticalGalaxyRe is Re in kpc
//OldStarFeedbackAlpha  alpha -19


  int i, j, k;
  float a=0;   

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1,
    TimeUnits = 1.0, VelocityUnits = 1.0;
  double MassUnits=1.0 ;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  float EnergyUnits;
  EnergyUnits = MassUnits* POW(VelocityUnits, 2.0);

  a=(EllipticalGalaxyRe*kpc_cm/1.8153)/LengthUnits;  // in code unit


  FLOAT r, x, y = 0, z = 0, rho_star=0;

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {

        /* Compute position */

        x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
        if (GridRank > 1)
          y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
        if (GridRank > 2)
          z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

          /* Find distance from center (0.5, 0.5, 0.5). */

        r = sqrt(POW(fabs(x-0.5), 2) +
                   POW(fabs(y-0.5), 2) +
                   POW(fabs(z-0.5), 2) );
        r = max(r, 0.1*CellWidth[0][0]);
        rho_star=ClusterSMBHBCG*SolarMass*1.0e11*a/(2.0*pi*r*POW((r+a)*LengthUnits,3)); //both r and a are in code units
        BaryonField[TENum][GRIDINDEX_NOGHOST(i,j,k)] += (0.000475*rho_star*(dtFixed*TimeUnits)*SNIaFeedbackEnergy  //SNIa
                                                        +rho_star*(dtFixed*TimeUnits)*OldStarFeedbackAlpha*1.0e-19*POW(300.0*1.0e5,2)) //Stellar
                                                        /(EnergyUnits*BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)]);   
        BaryonField[DensNum][GRIDINDEX_NOGHOST(i,j,k)] += rho_star*(dtFixed*TimeUnits)*OldStarFeedbackAlpha*1.0e-19/DensityUnits;
  }

  return SUCCESS;

}
 
