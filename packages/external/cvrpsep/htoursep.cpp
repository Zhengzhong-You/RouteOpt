/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include "memmod.h"
#include "basegrph.h"
#include "cnstrmgr.h"
#include "compress.h"
#include "capsep.h"
#include "newhtour.h"
#include "htoursep.h"

void HTOURSEP_SeparateHTours(int NoOfCustomers,
                             int *Demand,
                             int CAP,
                             int NoOfEdges,
                             int *EdgeTail,
                             int *EdgeHead,
                             double *EdgeX,
                             CnstrMgrPointer CMPExistingCuts,
                             int MaxNoOfCuts,
                             double *MaxViolation,
                             CnstrMgrPointer CutsCMP)
{
  int i,j,k;
  int NoOfV1Cuts;
  int ShrunkGraphCustNodes;
  int *SuperDemand;
  double *XInSuperNode;
  double **XMatrix;
  double **SMatrix;
  ReachPtr SupportPtr;
  ReachPtr V1CutsPtr;
  ReachPtr SAdjRPtr;
  ReachPtr SuperNodesRPtr;

  ReachInitMem(&SupportPtr,NoOfCustomers+1);
  ReachInitMem(&SAdjRPtr,NoOfCustomers+1);
  ReachInitMem(&SuperNodesRPtr,NoOfCustomers+1);

  SuperDemand = MemGetIV(NoOfCustomers+1);
  XInSuperNode = MemGetDV(NoOfCustomers+1);

  SMatrix = MemGetDM(NoOfCustomers+2,NoOfCustomers+2);
  XMatrix = MemGetDM(NoOfCustomers+2,NoOfCustomers+2);
  for (i=1; i<=NoOfCustomers+1; i++)
  for (j=1; j<=NoOfCustomers+1; j++)
  XMatrix[i][j] = 0.0;

  for (i=1; i<=NoOfEdges; i++)
  {
    ReachAddForwArc(SupportPtr,EdgeTail[i],EdgeHead[i]);
    ReachAddForwArc(SupportPtr,EdgeHead[i],EdgeTail[i]);

    XMatrix[EdgeTail[i]][EdgeHead[i]] = EdgeX[i];
    XMatrix[EdgeHead[i]][EdgeTail[i]] = EdgeX[i];
  }


  V1CutsPtr = NULL;
  CAPSEP_GetOneVehicleCapCuts(CMPExistingCuts,
                              &V1CutsPtr,
                              &NoOfV1Cuts);

  COMPRESS_ShrinkGraph(SupportPtr,
                       NoOfCustomers,
                       XMatrix,
                       SMatrix,
                       NoOfV1Cuts,
                       V1CutsPtr,
                       SAdjRPtr,
                       SuperNodesRPtr,
                       &ShrunkGraphCustNodes);

  ReachFreeMem(&V1CutsPtr);

  /* Compute data of supernodes */
  for (i=1; i<=ShrunkGraphCustNodes; i++)
  {
    XInSuperNode[i] = SMatrix[i][i];

    SuperDemand[i] = 0;
    for (j=1; j<=SuperNodesRPtr->LP[i].CFN; j++)
    {
      k = SuperNodesRPtr->LP[i].FAL[j];
      SuperDemand[i] += Demand[k];
    }
  }

  /* Ready to call Hypotour separation */

  NEWHTOUR_ComputeHTours(SupportPtr,
                         NoOfCustomers,
                         Demand,
                         CAP,
                         XMatrix,
                         SMatrix,
                         SuperNodesRPtr,
                         SAdjRPtr,
                         SuperDemand,
                         ShrunkGraphCustNodes,
                         CutsCMP,
                         MaxNoOfCuts,
                         MaxViolation);

  MemFree(SuperDemand);
  MemFree(XInSuperNode);

  MemFreeDM(SMatrix,NoOfCustomers+2);
  MemFreeDM(XMatrix,NoOfCustomers+2);

  ReachFreeMem(&SupportPtr);
  ReachFreeMem(&SAdjRPtr);
  ReachFreeMem(&SuperNodesRPtr);
}

