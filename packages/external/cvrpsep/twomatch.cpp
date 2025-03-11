/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "memmod.h"
#include "basegrph.h"
#include "cnstrmgr.h"
#include "sort.h"
#include "mxf.h"

void TWOMATCH_ComputeViolation(ReachPtr SupportPtr,
                               int NoOfCustomers,
                               double **XMatrix,
                               int *IntList,
                               int IntListSize,
                               int *ExtList,
                               int ExtListSize,
                               double *Violation)
{
  int i,j,k;
  int TSize;
  double TSum,OtherBoundary;
  char *InHandle;

  InHandle = MemGetCV(NoOfCustomers+2);
  for (i=1; i<=NoOfCustomers+1; i++) InHandle[i] = 0;

  for (i=1; i<=IntListSize; i++)
  {
    j = IntList[i];
    InHandle[j] = 1;
  }

  TSize = ExtListSize / 2;
  TSum = 0.0;
  OtherBoundary = 0.0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    if (InHandle[i] == 0) continue;

    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (InHandle[k] == 0)
      {
        OtherBoundary += XMatrix[i][k];
      }
    }
  }

  for (i=1; i<=TSize; i++)
  {
    j = ExtList[(2*i)-1];
    k = ExtList[2*i];
    TSum += XMatrix[j][k];
    OtherBoundary -= XMatrix[j][k];
  }

  *Violation = TSum - OtherBoundary - TSize + 1.0;

  MemFree(InHandle);
}

void TWOMATCH_GetCutNodeSet(ReachPtr RPtr,
                            int Source,
                            int *NodeList,
                            int *NodeListSize)
{
  int j,k;
  int Node,ReachedNr,ReachedListSize;

  ReachedListSize = 1;
  NodeList[1] = Source;

  ReachedNr = 0;
  while (ReachedNr < ReachedListSize)
  {
    ReachedNr++;
    Node = NodeList[ReachedNr];
    for (j=1; j<=RPtr->LP[Node].CFN; j++)
    {
      k = RPtr->LP[Node].FAL[j];
      NodeList[++ReachedListSize] = k;
    }
  }

  *NodeListSize = ReachedListSize;
}


void TWOMATCH_CheckConnectedHandle(ReachPtr SupportPtr,
                                   int NoOfCustomers,
                                   int *HandleList,
                                   int HandleListSize,
                                   char *Connected)
{
  int i,j,k,Node;
  int ReachedNr,ReachedListSize;
  char *InHandle, *Reached;
  int *ReachedList;

  InHandle = MemGetCV(NoOfCustomers+1);
  Reached  = MemGetCV(NoOfCustomers+1);

  ReachedList = MemGetIV(HandleListSize+1);

  for (i=1; i<=NoOfCustomers; i++) InHandle[i] = 0;
  for (i=1; i<=HandleListSize; i++)
  {
    j = HandleList[i];
    InHandle[j] = 1;
  }

  for (i=1; i<=NoOfCustomers; i++) Reached[i] = 0;

  j = HandleList[1];
  Reached[j] = 1;
  ReachedList[1] = j;
  ReachedListSize = 1;

  ReachedNr = 0;
  while (ReachedNr < ReachedListSize)
  {
    ReachedNr++;
    Node = ReachedList[ReachedNr];
    for (j=1; j<=SupportPtr->LP[Node].CFN; j++)
    {
      k = SupportPtr->LP[Node].FAL[j];
      if (k <= NoOfCustomers)
      if ((InHandle[k]) && (Reached[k] == 0))
      {
        ReachedList[++ReachedListSize] = k;
        Reached[k] = 1;
      }
    }
  }

  if (ReachedListSize == HandleListSize)
  *Connected = 1;
  else
  *Connected = 0;

  MemFree(InHandle);
  MemFree(Reached);
  MemFree(ReachedList);
}


void TWOMATCH_ExactTwoMatchings(ReachPtr SupportPtr,
                                int NoOfCustomers,
                                char *DepotEdgeBound,
                                double **XMatrix,
                                CnstrMgrPointer CutsCMP)
{
  char OddCut;
  char ConnectedHandle;
  char HandleOK;
  int i,j,k,Tail,Head;
  int NoOfNewNodes,NewNodeNr,FlowEdgesUB;
  int HandleSize;
  int TSize;
  int OriginalNodes,TotalNodes;
  double XVal,Scale,RHS;
  int ArcCap,MaxArcCap;
  int NodeListSize;
  int IntListSize, ExtListSize;
  char *OddNode, *OnHandleSide, *InTooth;
  int *NewNodeTail, *NewNodeHead;
  double *CutValue;
  int *NextOnPath;
  int *NodeList;
  int *IntList, *ExtList;
  MaxFlowPtr MXFPtr;
  ReachPtr RPtr;

  OriginalNodes = NoOfCustomers+1;

  NoOfNewNodes = 0;
  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (((k > i) && (k <= NoOfCustomers)) ||
         ((k == NoOfCustomers+1) && (DepotEdgeBound[i] == 1)))
      {
        NoOfNewNodes++;
      }
    }
  }

  /* (NewNodeTail[i],NewNodeHead[i]) is the original edge represented
     by new node #i) */

  NewNodeTail = MemGetIV(OriginalNodes+NoOfNewNodes+1);
  NewNodeHead = MemGetIV(OriginalNodes+NoOfNewNodes+1);

  NewNodeNr = 0;
  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (((k > i) && (k <= NoOfCustomers)) ||
         ((k == NoOfCustomers+1) && (DepotEdgeBound[i] == 1)))
      {
        NewNodeNr++;
        NewNodeTail[OriginalNodes+NewNodeNr] = i;
        NewNodeHead[OriginalNodes+NewNodeNr] = k;
      }
    }
  }

  TotalNodes = OriginalNodes + NoOfNewNodes;

  OddNode = MemGetCV(TotalNodes+1);
  OnHandleSide = MemGetCV(TotalNodes+1);
  InTooth = MemGetCV(TotalNodes+1);

  CutValue = MemGetDV(TotalNodes+1);
  NextOnPath = MemGetIV(TotalNodes+1);

  IntList = MemGetIV(TotalNodes+1);
  ExtList = MemGetIV(TotalNodes+1);


  /* (all new nodes are odd) */
  for (i=OriginalNodes+1; i<=TotalNodes; i++) OddNode[i] = 1;

  /* (compute odd/even nodes among the original customers) */
  for (i=1; i<=OriginalNodes; i++) OddNode[i] = 0;

  for (i=OriginalNodes+1; i<=TotalNodes; i++)
  {
    j = NewNodeTail[i];
    OddNode[j] = !OddNode[j];
  }

  FlowEdgesUB = (2 * NoOfNewNodes) + NoOfCustomers;
  FlowEdgesUB *= 2;
  MXF_InitMem(&MXFPtr,TotalNodes,FlowEdgesUB);
  MXF_ClearNodeList(MXFPtr);
  MXF_SetNodeListSize(MXFPtr,TotalNodes);
  MXF_ClearArcList(MXFPtr);

  Scale = 1000.0;
  MaxArcCap = 1000;

  for (i=OriginalNodes+1; i<=TotalNodes; i++)
  {
    Tail = NewNodeTail[i];
    Head = NewNodeHead[i];

    XVal = XMatrix[Tail][Head];
    ArcCap = (int)(XVal * Scale);
    if (ArcCap < 0) ArcCap = 0;
    if (ArcCap > MaxArcCap) ArcCap = MaxArcCap;

    if (MaxArcCap > ArcCap)
    {
      MXF_AddArc(MXFPtr,Tail,i,MaxArcCap-ArcCap);
      MXF_AddArc(MXFPtr,i,Tail,MaxArcCap-ArcCap);
    }

    if (ArcCap > 0)
    {
      MXF_AddArc(MXFPtr,Head,i,ArcCap);
      MXF_AddArc(MXFPtr,i,Head,ArcCap);
    }
  }

  /* (add edges to/from the depot) */
  /* (edges that have been split are already added in the above) */
  for (j=1; j<=SupportPtr->LP[NoOfCustomers+1].CFN; j++)
  {
    k = SupportPtr->LP[NoOfCustomers+1].FAL[j];

    /* (the edge for bound=1 has already been split and added) */
    if (DepotEdgeBound[k] == 1) continue;

    /* (add only the original capacity edge; 2-edges to/from the depot
        are not complemented) */

    XVal = XMatrix[NoOfCustomers+1][k];
    ArcCap = (int)(XVal * Scale);

    if (ArcCap > 0)
    {
      MXF_AddArc(MXFPtr,NoOfCustomers+1,k,ArcCap);
      MXF_AddArc(MXFPtr,k,NoOfCustomers+1,ArcCap);
    }
  }

  MXF_CreateMates(MXFPtr);

  MXF_ComputeGHCutTree(MXFPtr,
                       NoOfCustomers+1,
                       CutValue,
                       NextOnPath);

  ReachInitMem(&RPtr,TotalNodes);
  for (i=1; i<=TotalNodes; i++)
  {
    if (i != (NoOfCustomers+1))
    ReachAddForwArc(RPtr,NextOnPath[i],i);
  }

  NodeList = MemGetIV(TotalNodes+1);

  for (i=1; i<=TotalNodes; i++)
  {
    if (i == NoOfCustomers+1) continue;

    TWOMATCH_GetCutNodeSet(RPtr,i,NodeList,&NodeListSize);

    for (j=1; j<=TotalNodes; j++) OnHandleSide[j] = 0;
    for (j=1; j<=NodeListSize; j++) OnHandleSide[NodeList[j]] = 1;

    OddCut = 0;
    HandleSize = 0;

    for (j=1; j<=NodeListSize; j++)
    {
      if (NodeList[j] <= NoOfCustomers)
      HandleSize++;
      if (OddNode[NodeList[j]])
      OddCut = !OddCut;
    }

    if ((OddCut) && (HandleSize >= 3))
    {
      /* (setup the inequality) */
      IntListSize = 0;
      ExtListSize = 0;
      TSize = 0;

      for (j=1; j<=NodeListSize; j++)
      if (NodeList[j] <= NoOfCustomers)
      {
        IntList[++IntListSize] = NodeList[j];
      }

      if (IntListSize >= 3)
      {
        TWOMATCH_CheckConnectedHandle(SupportPtr,
                                      NoOfCustomers,
                                      IntList,
                                      IntListSize,
                                      &ConnectedHandle);
        HandleOK = ConnectedHandle;
      }
      else
      {
        HandleOK = 0;
      }

      if (HandleOK)
      {
        for (j=1; j<=OriginalNodes; j++) InTooth[j] = 0;

        for (j=OriginalNodes+1; j<=TotalNodes; j++)
        {
          Tail = NewNodeTail[j];
          Head = NewNodeHead[j];
          if ((OnHandleSide[Tail] != OnHandleSide[j]) &&
              (OnHandleSide[Tail] != OnHandleSide[Head]))
          {
            if ((InTooth[Tail] == 0) && (InTooth[Head] == 0))
            {
              TSize++;
              ExtList[++ExtListSize] = Tail;
              ExtList[++ExtListSize] = Head;

              /* (the depot may be included in more than one tooth) */
              if (Tail <= NoOfCustomers) InTooth[Tail] = 1;
              if (Head <= NoOfCustomers) InTooth[Head] = 1;
            }
          }
        }

      } /* (HandleOK) */

      if ((HandleOK) && (TSize >= 3) && ((TSize % 2) == 1))
      {
        RHS = IntListSize + ((TSize - 1) / 2);

        CMGR_AddExtCnstr(CutsCMP,
                         CMGR_CT_TWOMATCHING,
                         0,
                         IntListSize,
                         IntList,
                         ExtListSize,
                         ExtList,
                         RHS);

      }
    }
  }

  MemFree(NewNodeTail);
  MemFree(NewNodeHead);

  MemFree(OddNode);
  MemFree(OnHandleSide);
  MemFree(InTooth);

  MemFree(CutValue);
  MemFree(NextOnPath);

  MemFree(IntList);
  MemFree(ExtList);

  MemFree(NodeList);

  MXF_FreeMem(MXFPtr);
  ReachFreeMem(&RPtr);
}

