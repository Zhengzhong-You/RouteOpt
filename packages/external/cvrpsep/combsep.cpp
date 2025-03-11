/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include "memmod.h"
#include "basegrph.h"
#include "cnstrmgr.h"
#include "strcomb.h"

void COMBSEP_SeparateCombs(int NoOfCustomers,
                           int *Demand,
                           int CAP,
                           int QMin,
                           int NoOfEdges,
                           int *EdgeTail,
                           int *EdgeHead,
                           double *EdgeX,
                           int MaxNoOfCuts,
                           double *MaxViolation,
                           CnstrMgrPointer CutsCMP)
{
  int i,j;
  double **XMatrix;
  ReachPtr SupportPtr;

  ReachInitMem(&SupportPtr,NoOfCustomers+1);

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

  /* Ready to call comb separation */

  STRCOMB_MainStrengthenedCombs(SupportPtr,
                                NoOfCustomers,
                                CAP,
                                Demand,
                                QMin,
                                XMatrix,
                                MaxNoOfCuts,
                                MaxViolation,
                                CutsCMP);

  MemFreeDM(XMatrix,NoOfCustomers+2);
  ReachFreeMem(&SupportPtr);
}

