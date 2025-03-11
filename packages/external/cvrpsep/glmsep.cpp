/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include "memmod.h"
#include "basegrph.h"
#include "cnstrmgr.h"
#include "mxf.h"

void GLM_IdentifySingleSet(ReachPtr SupportPtr,
                           int *SmallGamma,
                           int BigGamma,
                           int NoOfCustomers,
                           double **XMatrix,
                           int *NodeList,
                           int *NodeListSize)
{
  int i,j,k,DepotIdx;
  int GraphNodes,ArcCap,CapFromSource;
  double XVal,XijFactor,ScaleFactor,WeightedSum,MaxFlowValue;
  double SGi, SGj, SGk, BG;
  int *SourceNodeList;
  double *DepotEdgeXVal;
  int SourceNodeListSize;
  MaxFlowPtr MXFPtr;

  ScaleFactor = 1000.0;
  BG = BigGamma;

  GraphNodes = (NoOfCustomers + 2);
  MXF_InitMem(&MXFPtr,GraphNodes,GraphNodes*5);
  MXF_ClearNodeList(MXFPtr);
  MXF_SetNodeListSize(MXFPtr,GraphNodes);
  MXF_ClearArcList(MXFPtr);

  SourceNodeList = MemGetIV(NoOfCustomers+2);
  DepotEdgeXVal = MemGetDV(NoOfCustomers+1);

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if ((k <= NoOfCustomers) && (k > i))
      {
        XVal = XMatrix[i][k];

        SGi = SmallGamma[i];
        SGk = SmallGamma[k];

        XijFactor = 1.0 - ((SGi + SGk) / BG);
        XVal *= XijFactor;

        ArcCap = (int)(XVal * ScaleFactor);

        MXF_AddArc(MXFPtr,i,k,ArcCap);
        MXF_AddArc(MXFPtr,k,i,ArcCap);
      }
    }
  }

  DepotIdx = NoOfCustomers + 1;

  for (k=1; k<=NoOfCustomers; k++) DepotEdgeXVal[k] = 0.0;

  for (j=1; j<=SupportPtr->LP[DepotIdx].CFN; j++)
  {
    k = SupportPtr->LP[DepotIdx].FAL[j];
    DepotEdgeXVal[k] = XMatrix[DepotIdx][k];
  }

  CapFromSource = 0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    SGi = SmallGamma[i];
    XVal = (1.0 - (SGi / BG)) * DepotEdgeXVal[i];

    WeightedSum = 0.0;
    for (k=1; k<=SupportPtr->LP[i].CFN; k++)
    {
      j = SupportPtr->LP[i].FAL[k];
      if (j <= NoOfCustomers)
      {
        SGj = SmallGamma[j];
        WeightedSum += (SGj / BG) * XMatrix[i][j];
      }
    }

    XVal -= WeightedSum;

    ArcCap = (int)(XVal * ScaleFactor);

    if (ArcCap > 0)
    {
      MXF_AddArc(MXFPtr,i,NoOfCustomers+2,ArcCap);
    }
    else
    {
      ArcCap = -ArcCap;
      ArcCap++;
      MXF_AddArc(MXFPtr,NoOfCustomers+1,i,ArcCap);
      CapFromSource += ArcCap;
    }
  }

  if (CapFromSource > 0)
  {
    MXF_CreateMates(MXFPtr);
    MXF_SolveMaxFlow(MXFPtr,1,NoOfCustomers+1,NoOfCustomers+2,&MaxFlowValue,0,
                      &SourceNodeListSize,SourceNodeList);

    SourceNodeListSize--;
  }
  else
  {
    SourceNodeListSize = 0;
  }

  for (i=1; i<=SourceNodeListSize; i++) NodeList[i] = SourceNodeList[i];
  *NodeListSize = SourceNodeListSize;

  MemFree(SourceNodeList);
  MemFree(DepotEdgeXVal);

  MXF_FreeMem(MXFPtr);
}

void GLMSEP_SeparateGLM(int NoOfCustomers,
                        int *Demand,
                        int CAP,
                        int NoOfEdges,
                        int *EdgeTail,
                        int *EdgeHead,
                        double *EdgeX,
                        int *CustList,
                        int *CustListSize,
                        double *Violation)
{
  int i,j;
  double LHS,RHS;
  char *InNodeSet;
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

  /* Ready to call GLM separation */
  GLM_IdentifySingleSet(SupportPtr,Demand,CAP,NoOfCustomers,
                        XMatrix,CustList,CustListSize);

  if (*CustListSize > 0)
  {
    InNodeSet = MemGetCV(NoOfCustomers+1);
    for (i=1; i<=NoOfCustomers; i++) InNodeSet[i] = 0;
    for (i=1; i<=*CustListSize; i++)
    {
      j = CustList[i];
      InNodeSet[j] = 1;
    }

    LHS = 0.0;
    RHS = 0.0;

    for (i=1; i<=NoOfCustomers; i++)
    if (InNodeSet[i])
    {
      LHS += (CAP * XMatrix[i][NoOfCustomers+1]);
      RHS += (2 * Demand[i]);
    }

    for (i=1; i<=NoOfCustomers; i++)
    if (InNodeSet[i])
    {
      for (j=1; j<=NoOfCustomers; j++)
      if (InNodeSet[j] == 0)
      {
        LHS += ((CAP - (2*Demand[j])) * XMatrix[i][j]);
      }
    }

    *Violation = (RHS - LHS) / (1.0 * CAP);

    if (*Violation < 0.02)
    {
      *CustListSize = 0;
      *Violation = 0.0;
    }

    MemFree(InNodeSet);
  }

  MemFreeDM(XMatrix,NoOfCustomers+2);
  ReachFreeMem(&SupportPtr);
}

