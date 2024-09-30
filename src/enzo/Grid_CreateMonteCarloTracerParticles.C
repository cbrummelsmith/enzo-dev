#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h" 

void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);

int grid::CreateMonteCarloTracerParticles()
{

  int MCTP_CREATION_MODE = 0;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int groupID = 0;
  int size = 1;
  int i, j, k, n, index;
  FLOAT pos[3];
  MonteCarloTracerParticle *newMC;

  printf("\n%s\n", "CreateMonteCarloTracerParticles...");

  this->AllocateMonteCarloTracerParticleData();

  printf("\n%s\n", "AllocatedMonteCarloTracerParticles...");  

  /* Create N particles per cel */

  if (MCTP_CREATION_MODE == 0) {

    // Loop over all active cells 
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      pos[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][0];
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        pos[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
          pos[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][0]; 

          // /* DEBUG */
          // if (i != GridStartIndex[0] || j != GridStartIndex[1] || k != GridStartIndex[2])
          //   continue;

          // /* DEBUG */
          // if (i != GridStartIndex[0] || j != GridEndIndex[1] || k != GridStartIndex[2])
          //   continue;        

          // /* DEBUG */
          // if (j != GridStartIndex[1])
          //   continue;          

          // if (i != GridEndIndex[0] || j != GridEndIndex[1] || k != GridEndIndex[2])
          //   continue;            

          // if (i % 2 != 0 || j % 2 != 0 || k % 2 != 0)
          //   continue;

          /* END DEBUG */

          index = (k * GridDimension[1] + j) * GridDimension[0] + i;

          // Add MC tracers to this cell 
          for (n = 0; n < NumberOfMonteCarloTracerParticlesPerCell; n++) {
            //printf("\n(i,j,k,n) = (%d,%d,%d,%d) New Particle", i, j, k, n);  
            newMC = new MonteCarloTracerParticle(this, index, groupID, GridLevel, Time, pos);
            //printf("\nnewMC %p", newMC);            
            InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[index], newMC);
            //printf("\nInserted %p", newMC);            
            TotalNumberOfMonteCarloTracerParticles++;
          }
        }
      }
    }
  } /* END MCTP_CREATION_MODE == 0 */

  /* Create particles per cell based on gas density. 
     N particles in max density cells. 
     n/cell scales down with density */

  else if (MCTP_CREATION_MODE == 1) {

    /* compute the field size */
  
    int size = 1;
    float min_density = huge_number;
    int* MCTP_per_cell;

    for (int dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
    }

    MCTP_per_cell = new int[size];

    /* Find fields: density, total energy, velocity1-3. */
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                        Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

    /* Find max density */
    for (i = 0; i < size; i++)
        if (BaryonField[DensNum][i] < min_density)
          min_density = BaryonField[DensNum][i];

    /* Compute number of particles for each cell based on density ratio*/
    for (i = 0; i < size; i++)
      MCTP_per_cell[i] = int(BaryonField[DensNum][i]/min_density * NumberOfMonteCarloTracerParticlesPerCell);

    // Loop over all active cells 
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      pos[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][0];
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        pos[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
          pos[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][0]; 

          index = (k * GridDimension[1] + j) * GridDimension[0] + i;

          // Add MC tracers to this cell 
          for (n = 0; n < MCTP_per_cell[index]; n++) {
            newMC = new MonteCarloTracerParticle(this, index, groupID, GridLevel, Time, pos);
            InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[index], newMC);
            TotalNumberOfMonteCarloTracerParticles++;
          }

        } // End i
      } // End j
    } // End k
  } /* END MCTP_CREATION_MODE == 1*/

  printf("\nTotalNumberOfMonteCarloTracerParticles %i\n", TotalNumberOfMonteCarloTracerParticles);
  printf("%s\n", "CreateMonteCarloTracerParticles Done.");

  return SUCCESS;
}
