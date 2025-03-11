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
#include "strngcmp.h"
#include "sort.h"
#include "blocks.h"
#include "cutbase.h"
#include "cnstrmgr.h"
#include "grsearch.h"
#include "twomatch.h"

void STRCOMB_ComputeBoundaryLHS(ReachPtr SupportPtr,
                                int NoOfCustomers,
                                double **XMatrix,
                                int NoOfTeeth,
                                int *IntList,
                                int IntListSize,
                                int *ExtList,
                                int ExtListSize,
                                double *LHS)
{
  int i,j,k,t;
  int MinIdx,MaxIdx;
  char **InTooth;

  *LHS = 0.0;

  InTooth = MemGetCM(NoOfCustomers+2,NoOfTeeth+1);
  for (i=1; i<=NoOfCustomers+1; i++)
  for (t=0; t<=NoOfTeeth; t++) /* Handle = tooth nr. 0 */
  InTooth[i][t] = 0;

  for (i=1; i<=IntListSize; i++)
  {
    j = IntList[i];
    InTooth[j][0] = 1;
  }

  for (t=1; t<=NoOfTeeth; t++)
  {
    MinIdx = ExtList[t];

    if (t == NoOfTeeth)
    MaxIdx = ExtListSize;
    else
    MaxIdx = ExtList[t+1] - 1;

    for (i=MinIdx; i<=MaxIdx; i++)
    {
      j = ExtList[i];
      InTooth[j][t] = 1;
    }
  }

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (k > i)
      {
        for (t=0; t<=NoOfTeeth; t++)
        {
          if ((InTooth[i][t]) != (InTooth[k][t]))
          (*LHS) += XMatrix[i][k];
        }
      }
    }
  }

  MemFreeCM(InTooth,NoOfCustomers+2);
}

void STRCOMB_ComputeStrCombRHS(int NoOfCustomers,
                               int NoOfTeeth,
                               int *Demand,
                               int CAP,
                               int *IntList,
                               int IntListSize,
                               int *ExtList,
                               int ExtListSize,
                               int *RHS)
{
  char DepotInTooth;
  int i,j,t;
  int MinIdx,MaxIdx;
  int TotalDemand;
  int *Q, *QCalc, *QCap, *QK;
  char **InTooth;

  Q     = MemGetIV(4);
  QCalc = MemGetIV(4);
  QCap  = MemGetIV(4);
  QK    = MemGetIV(4);

  InTooth = MemGetCM(NoOfCustomers+2,NoOfTeeth+1);

  for (i=1; i<=NoOfCustomers+1; i++)
  for (t=0; t<=NoOfTeeth; t++) /* Handle = tooth nr. 0 */
  InTooth[i][t] = 0;

  for (i=1; i<=IntListSize; i++)
  {
    j = IntList[i];
    InTooth[j][0] = 1;
  }

  for (t=1; t<=NoOfTeeth; t++)
  {
    MinIdx = ExtList[t];

    if (t == NoOfTeeth)
    MaxIdx = ExtListSize;
    else
    MaxIdx = ExtList[t+1] - 1;

    for (i=MinIdx; i<=MaxIdx; i++)
    {
      j = ExtList[i];
      InTooth[j][t] = 1;
    }
  }

  TotalDemand = 0;
  for (i=1; i<=NoOfCustomers; i++) TotalDemand += Demand[i];

  *RHS = 0;

  for (t=1; t<=NoOfTeeth; t++)
  {
    DepotInTooth = InTooth[NoOfCustomers+1][t];

    for (i=1; i<=3; i++) Q[i] = 0;

    for (i=1; i<=NoOfCustomers; i++)
    {
      if (InTooth[i][t])
      { /* Node i is in the tooth */
        if (InTooth[i][0])
        {
          /* Node i is in the handle */
          Q[1] += Demand[i];
        }
        else
        {
          /* Node i is not in the handle */
          Q[2] += Demand[i];
        }

        Q[3] += Demand[i];
      }
    }

    for (i=1; i<=3; i++)
    QCalc[i] = Q[i];

    if (DepotInTooth)
    {
      QCalc[2] = TotalDemand - Q[2];
      QCalc[3] = TotalDemand - Q[3];
    }

    for (i=1; i<=3; i++)
    {
      QCap[i] = CAP;
      QK[i] = 1;

      while (QCap[i] < QCalc[i])
      {
        QCap[i] += CAP;
        (QK[i])++;
      }
    }

    (*RHS) += (QK[1] + QK[2] + QK[3]);
  }

  MemFree(Q);
  MemFree(QCalc);
  MemFree(QCap);
  MemFree(QK);

  MemFreeCM(InTooth,NoOfCustomers+2);
}

void STRCOMB_ExpandTooth(ReachPtr SupportPtr,
                         int NoOfCustomers,
                         int NoOfTeeth,
                         int ThisToothNr,
                         int *Demand,
                         int CAP,
                         double *NodeBoundary,
                         char *InHandle,
                         char **InTooth,
                         double **XMatrix,
                         ReachPtr ToothRPtr,
                         int *ParamBestList,
                         int *ParamBestListSize,
                         double *LHS,
                         int *RHS)
{
  char DepotInTooth,ExpandTooth;
  char TeethIntersect,IntersectionInHandle;
  int i,j,k,t;
  int TotalDemand;
  int NodeListSize,BestNodeListSize,InitNodeListSize,BestNode,BestXNode;
  int Stage,BestStage = 0;
  double XVal,BestX,ToothDelta,Slack,BestSlack;
  double NewToothDelta,NewSlack,NewBestSlack;
  char *Reached, *Selectable, *InSet;
  char *TouchedTooth;
  int *Q, *QCalc, *QCap, *QK;
  int *NewQ, *NewQCalc, *NewQCap, *NewQK;
  int *NodeStage, *NodeList;
  double *XValue;

  Q     = MemGetIV(4);
  QCalc = MemGetIV(4);
  QCap  = MemGetIV(4);
  QK    = MemGetIV(4);

  NewQ     = MemGetIV(4);
  NewQCalc = MemGetIV(4);
  NewQCap  = MemGetIV(4);
  NewQK    = MemGetIV(4);

  TouchedTooth = MemGetCV(NoOfTeeth+1);

  NodeList = MemGetIV(NoOfCustomers+2);
  NodeStage = MemGetIV(NoOfCustomers+2);
  InSet = MemGetCV(NoOfCustomers+2);
  Reached = MemGetCV(NoOfCustomers+1);
  Selectable = MemGetCV(NoOfCustomers+2);
  XValue = MemGetDV(NoOfCustomers+1);

  Selectable[NoOfCustomers+1] = 0;
  for (i=1; i<=NoOfCustomers; i++) Selectable[i] = 1;

  Stage=ToothRPtr->LP[ThisToothNr].CFN; /* Initial size of the tooth */
  for (i=1; i<=NoOfCustomers+1; i++) NodeStage[i] = NoOfCustomers+2;

  for (i=1; i<=NoOfCustomers+1; i++) InSet[i] = InTooth[i][ThisToothNr];

  for (j=1; j<=ToothRPtr->LP[ThisToothNr].CFN; j++)
  {
    k = ToothRPtr->LP[ThisToothNr].FAL[j];
    NodeStage[k] = Stage;
    Selectable[k] = 0;
  }

  for (t=1; t<=NoOfTeeth; t++) TouchedTooth[t] = 0;
  TouchedTooth[ThisToothNr] = 1;

  for (t=1; t<=NoOfTeeth; t++)
  {
    if (t == ThisToothNr) continue;

    TeethIntersect=0;
    for (j=1; j<=ToothRPtr->LP[t].CFN; j++)
    {
      k = ToothRPtr->LP[t].FAL[j];

      if (NodeStage[k] <= Stage) /* Intersection with current tooth */
      {
        TeethIntersect=1;
        IntersectionInHandle = InHandle[k];
        TouchedTooth[t] = 1;
        break;
      }
    }

    if (TeethIntersect)
    {
      for (j=1; j<=ToothRPtr->LP[t].CFN; j++)
      {
        k = ToothRPtr->LP[t].FAL[j];

        if (InHandle[k] != IntersectionInHandle)
        Selectable[k] = 0;
      }
    }

    /* If we don't permit intersection inside the handle */

    for (j=1; j<=ToothRPtr->LP[t].CFN; j++)
    {
      k = ToothRPtr->LP[t].FAL[j];
      if (InHandle[k])
      Selectable[k] = 0;
    }
  }

  /* OK: Now we don't get teeth that intersect both inside and */
  /* outside the handle. */

  for (i=1; i<=NoOfCustomers; i++) Reached[i] = 0;
  for (i=1; i<=NoOfCustomers; i++) XValue[i] = 0.0;

  DepotInTooth = InSet[NoOfCustomers+1];
  ToothDelta = 0.0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (k > i)
      {
        XVal = XMatrix[i][k];

        if ((InSet[i]) && (k <= NoOfCustomers))
        {
          XValue[k] += XVal;
          Reached[k] = 1;
        }

        if (InSet[k])
        XValue[i] += XVal;

        if ((InSet[k]) && (k <= NoOfCustomers))
        Reached[i] = 1;

        if (InSet[i] != InSet[k])
        {
          ToothDelta += XVal;
        }
      }
    }
  }

  for (i=1; i<=3; i++) Q[i] = 0;

  TotalDemand = 0;
  for (i=1; i<=NoOfCustomers; i++)
  {
    TotalDemand += Demand[i];

    if (InSet[i])
    { /* Node i is in the tooth */
      if (InHandle[i])
      {
        /* Node i is in the handle */
        Q[1] += Demand[i];
      }
      else
      {
        /* Node i is not in the handle */
        Q[2] += Demand[i];
      }

      Q[3] += Demand[i];
    }
  }

  for (i=1; i<=3; i++)
  QCalc[i] = Q[i];

  if (DepotInTooth)
  { /* The depot is outside the handle */
    QCalc[2] = TotalDemand - Q[2];
    QCalc[3] = TotalDemand - Q[3];
  }

  for (i=1; i<=3; i++)
  {
    QCap[i] = CAP;
    QK[i] = 1;

    while (QCap[i] < QCalc[i])
    {
      QCap[i] += CAP;
      (QK[i])++;
    }
  }

  *LHS = ToothDelta;
  *RHS = QK[1] + QK[2] + QK[3];
  Slack = *LHS - *RHS;

  NodeListSize = 0;
  for (i=1; i<=NoOfCustomers+1; i++)
  if (InSet[i])
  {
    NodeList[++NodeListSize] = i;
  }

  InitNodeListSize = NodeListSize;
  BestNodeListSize = NodeListSize;
  BestSlack = Slack;

  ExpandTooth = 1;
  while (ExpandTooth)
  {
    /* Add the customer with max x-value */

    BestX = -1.0;
    BestNode = 0;
    BestXNode = 0;

    NewBestSlack = 2.0 * (NoOfCustomers + 2);

    for (i=1; i<=NoOfCustomers; i++)
    {
      if (Selectable[i] == 0) continue;
      if (Reached[i] == 0) continue;

      if (XValue[i] > BestX)
      {
        BestX = XValue[i];
        BestXNode = i;
      }

      /* Calculate best slack */
      for (j=1; j<=3; j++) NewQ[j] = Q[j];

      if (InHandle[i])
      NewQ[1] += Demand[i];
      else
      NewQ[2] += Demand[i];

      NewQ[3] += Demand[i];

      for (j=1; j<=3; j++)
      {
        NewQCalc[j] = NewQ[j];
        NewQCap[j] = QCap[j];
        NewQK[j] = QK[j];
      }

      if (DepotInTooth)
      {
        NewQCalc[2] = TotalDemand - NewQ[2];
        NewQCalc[3] = TotalDemand - NewQ[3];
      }

      for (j=1; j<=3; j++)
      {
        while (NewQCap[j] < NewQCalc[j])
        {
          NewQCap[j] += CAP;
          (NewQK[j])++;
        }

        while ((NewQCap[j] - CAP) >= NewQCalc[j])
        {
          NewQCap[j] -= CAP;
          (NewQK[j])--;
        }
      }

      NewToothDelta = ToothDelta + (NodeBoundary[i] - (2.0 * XValue[i]));
      NewSlack = NewToothDelta - (NewQK[1] + NewQK[2] + NewQK[3]);

      if (NewSlack < (NewBestSlack - 0.001))
      {
        NewBestSlack = NewSlack;
        BestNode = i;
      }
    }

    if (BestNode == 0)
    {
      BestNode = BestXNode;
    }

    if (BestNode > 0)
    {
      BestX = XValue[BestNode];

      /* Include BestNode */
      InSet[BestNode] = 1;
      Selectable[BestNode] = 0;
      NodeList[++NodeListSize] = BestNode;
      NodeStage[BestNode] = ++BestStage;

      if (InHandle[BestNode])
      {
        /* Node BestNode is in the handle */
        Q[1] += Demand[BestNode];
      }
      else
      {
        /* Node BestNode is not in the handle */
        Q[2] += Demand[BestNode];
      }

      Q[3] += Demand[BestNode];

      for (i=1; i<=3; i++) QCalc[i] = Q[i];

      if (DepotInTooth)
      {
        QCalc[2] = TotalDemand - Q[2];
        QCalc[3] = TotalDemand - Q[3];
      }

      for (i=1; i<=3; i++)
      {
        while (QCap[i] < QCalc[i])
        {
          QCap[i] += CAP;
          (QK[i])++;
        }

        while ((QCap[i] - CAP) >= QCalc[i])
        {
          QCap[i] -= CAP;
          (QK[i])--;
        }
      }

      ToothDelta += (NodeBoundary[BestNode] - (2.0 * BestX));
      Slack = ToothDelta - (QK[1] + QK[2] + QK[3]);

      if ((Slack < (BestSlack - 0.001)) || ((*RHS % 2) == 0))
      if ((QK[1] + QK[2] + QK[3]) % 2 == 1) /* Odd sum */
      {
        BestNodeListSize = NodeListSize;
        BestSlack = Slack;

        *LHS = ToothDelta;
        *RHS = QK[1] + QK[2] + QK[3];
      }

      for (j=1; j<=SupportPtr->LP[BestNode].CFN; j++)
      {
        k = SupportPtr->LP[BestNode].FAL[j];
        if (k <= NoOfCustomers)
        {
          XVal = XMatrix[BestNode][k];
          XValue[k] += XVal;
          Reached[k] = 1;
        }
      }

      /* Identify new non-selectable nodes */
      for (t=1; t<=NoOfTeeth; t++)
      {
        if (t == ThisToothNr) continue;
        if (InTooth[BestNode][t] == 0) continue;

        if (TouchedTooth[t] == 1) continue;
        /* We already have a node from tooth t in ThisTooth */

        IntersectionInHandle = InHandle[BestNode];

        for (j=1; j<=ToothRPtr->LP[t].CFN; j++)
        {
          k = ToothRPtr->LP[t].FAL[j];

          if (InHandle[k] != IntersectionInHandle)
          {
            Selectable[k] = 0;
          }
        }

        TouchedTooth[t] = 1;
      }

    }

    if (BestNode == 0) break;
  }

  *ParamBestListSize = BestNodeListSize;
  for (i=1; i<=BestNodeListSize; i++)
  ParamBestList[i] = NodeList[i];

  MemFree(Reached);
  MemFree(Selectable);
  MemFree(InSet);
  MemFree(TouchedTooth);

  MemFree(Q);
  MemFree(QCalc);
  MemFree(QCap);
  MemFree(QK);

  MemFree(NewQ);
  MemFree(NewQCalc);
  MemFree(NewQCap);
  MemFree(NewQK);

  MemFree(NodeStage);
  MemFree(NodeList);

  MemFree(XValue);
}


void STRCOMB_ExpandToothTwoWays(ReachPtr SupportPtr,
                                int NoOfCustomers,
                                int NoOfTeeth,
                                int ThisToothNr,
                                int *Demand,
                                int CAP,
                                double *NodeBoundary,
                                char *InHandle,
                                char **InTooth,
                                double **XMatrix,
                                ReachPtr ToothRPtr,
                                double *LHS,
                                int *RHS)
{
  char DepotInTooth,AddDepot;
  int LocalRHS;
  double LocalLHS;
  int i,j,k,t;
  int ListSize, BestListSize;
  double Slack,BestSlack;
  int *List, *BestList;

  List = MemGetIV(NoOfCustomers+2);
  BestList = MemGetIV(NoOfCustomers+2);

  DepotInTooth = InTooth[NoOfCustomers+1][ThisToothNr];

  STRCOMB_ExpandTooth(SupportPtr,
                      NoOfCustomers,
                      NoOfTeeth,
                      ThisToothNr,
                      Demand,
                      CAP,
                      NodeBoundary,
                      InHandle,
                      InTooth,
                      XMatrix,
                      ToothRPtr,
                      BestList,
                      &BestListSize,
                      &LocalLHS,
                      &LocalRHS);

  *LHS = LocalLHS;
  *RHS = LocalRHS;
  BestSlack = LocalLHS - LocalRHS;

  if (DepotInTooth == 0)
  {
    /* Try with the depot inside the tooth */

    /* Check if it is valid to add the depot to this tooth */
    AddDepot = 1;

    for (t=1; t<=NoOfTeeth; t++)
    {
      if (t == ThisToothNr) continue;
      if (InTooth[NoOfCustomers+1][t] == 0) continue;

      /* The depot is in tooth t */

      /* if tooth t intersects with ThisToothNr inside the handle, */
      /* then we cannot add the depot to ThisToothNr outside the handle. */

      for (j=1; j<=ToothRPtr->LP[t].CFN; j++)
      {
        k = ToothRPtr->LP[t].FAL[j];

        if (InTooth[k][ThisToothNr])
        {
          AddDepot = 0;
          break;
        }
      }
    }

    if (AddDepot)
    {
      InTooth[NoOfCustomers+1][ThisToothNr] = 1;

      ListSize = 0;
      for (i=1; i<=NoOfCustomers+1; i++)
      if (InTooth[i][ThisToothNr]) List[++ListSize] = i;

      ReachSetForwList(ToothRPtr,List,ThisToothNr,ListSize);
      /* Now including the depot */

      STRCOMB_ExpandTooth(SupportPtr,
                          NoOfCustomers,
                          NoOfTeeth,
                          ThisToothNr,
                          Demand,
                          CAP,
                          NodeBoundary,
                          InHandle,
                          InTooth,
                          XMatrix,
                          ToothRPtr,
                          List,
                          &ListSize,
                          &LocalLHS,
                          &LocalRHS);

      Slack = LocalLHS - LocalRHS;

      if ((LocalRHS % 2) == 1)
      {
        if ((Slack < (BestSlack - 0.01)) || ((*RHS % 2) == 0))
        {
          BestSlack = Slack;

          *LHS = LocalLHS;
          *RHS = LocalRHS;

          for (i=1; i<=ListSize; i++) BestList[i] = List[i];
          BestListSize = ListSize;
        }
      }
    } /* AddDepot */
  }

  for (i=1; i<=NoOfCustomers+1; i++) InTooth[i][ThisToothNr] = 0;

  for (i=1; i<=BestListSize; i++)
  {
    j = BestList[i];
    InTooth[j][ThisToothNr] = 1;
  }

  ReachSetForwList(ToothRPtr,BestList,ThisToothNr,BestListSize);

  MemFree(List);
  MemFree(BestList);
}

void STRCOMB_GenerateHandlesV2(ReachPtr SupportPtr,
                               int NoOfCustomers,
                               double **XMatrix,
                               ReachPtr HandlesRPtr,
                               int MaxNoOfHandles,
                               int *FoundNoOfHandles)
{
  char ListFound;
  char TailInBlock,HeadInBlock;
  int i,j,k,Idx,EdgeNr,Tail,Head,Label,NodeSum;
  int TailComp,HeadComp;
  int BlockNr;
  int NoOfEdges,NoOfBlocks,NoOfHandles,NodeListSize;
  double XVal;
  int *EdgeHead, *EdgeTail, *NodeList;
  int *Index;
  int *NodeLabel;
  double *Value;
  int *CompNr;
  ReachPtr EdgesRPtr, BlocksRPtr;

  NoOfHandles=0;

  NoOfEdges = 0;
  for (i=1; i<NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if ((k > NoOfCustomers) || (k < i)) continue;

      NoOfEdges++;
    }
  }

  CompNr  = MemGetIV(NoOfCustomers+1);

  NodeList = MemGetIV(NoOfCustomers+1);
  NodeLabel = MemGetIV(NoOfCustomers+1);

  EdgeHead = MemGetIV(NoOfEdges+1);
  EdgeTail = MemGetIV(NoOfEdges+1);
  Index    = MemGetIV(NoOfEdges+1);
  Value    = MemGetDV(NoOfEdges+1);

  EdgeNr = 0;
  for (i=1; i<NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if ((k > NoOfCustomers) || (k < i)) continue;

      EdgeNr++;
      EdgeTail[EdgeNr] = i;
      EdgeHead[EdgeNr] = k;
      Index[EdgeNr]    = EdgeNr;

      XVal = XMatrix[i][k] - 0.5;
      if (XVal < 0.0) XVal = -XVal;
      Value[EdgeNr] = XVal;
    }
  }

  SortIndexDVInc(Index,Value,NoOfEdges);

  ReachInitMem(&EdgesRPtr,NoOfCustomers);
  ReachInitMem(&BlocksRPtr,NoOfEdges);

  for (i=1; i<=NoOfCustomers; i++) CompNr[i] = i;

  /* Components: */

  for (Idx=1; Idx<=NoOfEdges; Idx++)
  {
    if (NoOfHandles >= MaxNoOfHandles)
    {
      break;
    }

    EdgeNr = Index[Idx];

    Tail = EdgeTail[EdgeNr];
    Head = EdgeHead[EdgeNr];

    if (CompNr[Tail] != CompNr[Head])
    {
      /* Merge the two components into one */
      /* Move tail to head */

      TailComp = CompNr[Tail];
      HeadComp = CompNr[Head];

      NodeListSize = 0;
      for (i=1; i<=NoOfCustomers; i++)
      {
        if ((CompNr[i] == TailComp) || (CompNr[i] == HeadComp))
        {
          CompNr[i] = HeadComp;
          NodeList[++NodeListSize] = i;
        }
      }

      if (NodeListSize >= 3)
      {
        NoOfHandles++;
        GRSEARCH_AddSet(HandlesRPtr,
                        NoOfHandles,
                        NodeListSize,
                        NodeList,
                        0);
      }
    }
  }

  /* Blocks: */

  Label = 0;
  for (i=1; i<=NoOfCustomers; i++) NodeLabel[i] = 0;

  for (Idx=1; Idx<=NoOfEdges; Idx++)
  {
    if (NoOfHandles >= MaxNoOfHandles)
    {
      break;
    }

    Label++;
    EdgeNr = Index[Idx];

    Tail = EdgeTail[EdgeNr];
    Head = EdgeHead[EdgeNr];

    ReachAddForwArc(EdgesRPtr,Tail,Head);
    ReachAddForwArc(EdgesRPtr,Head,Tail);

    ComputeBlocks(EdgesRPtr,BlocksRPtr,NoOfCustomers,&NoOfBlocks);

    for (BlockNr=1; BlockNr<=NoOfBlocks; BlockNr++)
    {
      if (BlocksRPtr->LP[BlockNr].CFN <= 2) continue;

      TailInBlock = 0;
      HeadInBlock = 0;

      for (j=1; j<=BlocksRPtr->LP[BlockNr].CFN; j++)
      {
        k = BlocksRPtr->LP[BlockNr].FAL[j];
        if (k == Tail) TailInBlock = 1;
        if (k == Head) HeadInBlock = 1;
      }

      if ((TailInBlock) && (HeadInBlock))
      {
        NodeSum = 0;
        NodeListSize = 0;

        for (j=1; j<=BlocksRPtr->LP[BlockNr].CFN; j++)
        {
          k = BlocksRPtr->LP[BlockNr].FAL[j];
          NodeList[++NodeListSize] = k;
          NodeLabel[k] = Label;
          NodeSum += k;
        }

        GRSEARCH_CheckForExistingSet(HandlesRPtr,
                                     NoOfHandles,
                                     NodeLabel,
                                     Label,
                                     NodeSum,
                                     NodeListSize,
                                     &ListFound);

        if (ListFound == 0)
        {
          NoOfHandles++;
          GRSEARCH_AddSet(HandlesRPtr,
                          NoOfHandles,
                          NodeListSize,
                          NodeList,
                          0);
        }

        break;
      } /* Edge in block */
    }
  }

  *FoundNoOfHandles = NoOfHandles;

  MemFree(CompNr);

  MemFree(NodeList);
  MemFree(NodeLabel);

  MemFree(EdgeHead);
  MemFree(EdgeTail);
  MemFree(Index);
  MemFree(Value);

  ReachFreeMem(&EdgesRPtr);
  ReachFreeMem(&BlocksRPtr);
}

void STRCOMB_Shrink(ReachPtr SupportPtr,
                    int NoOfCustomers,
                    int *Demand,
                    int QMin,
                    char UseDepotMatch,
                    double **XMatrix,
                    int *CompNr, /* Node j is in component nr. CompNr[j]. */
                    int *CompDemand,
                    double *CompBoundary,
                    double **SMatrix,
                    ReachPtr ResultCompsRPtr,
                    int *ResultNoOfComponents)
{
  char ShrinkAgain;
  int i,j,k,Idx;
  int NoOfComponents,DepotCompNr;
  int CompI,CompJ,CompK;
  double XVal,EdgeSum;

  char *CVWrk1;
  int *IVWrk1, *IVWrk2, *IVWrk3, *IVWrk4;

  ReachPtr CmprsEdgesRPtr;
  ReachPtr CompsRPtr;

  CVWrk1 = MemGetCV(NoOfCustomers+2);
  IVWrk1 = MemGetIV(NoOfCustomers+2);
  IVWrk2 = MemGetIV(NoOfCustomers+2);
  IVWrk3 = MemGetIV(NoOfCustomers+2);
  IVWrk4 = MemGetIV(NoOfCustomers+2);

  ReachInitMem(&CmprsEdgesRPtr,NoOfCustomers+1);
  ReachInitMem(&CompsRPtr,NoOfCustomers+1);

  do
  {
    ShrinkAgain = 0;

    ReachClearForwLists(CompsRPtr);
    ComputeStrongComponents(CmprsEdgesRPtr,CompsRPtr,
                            &NoOfComponents,NoOfCustomers+1,
                            CVWrk1,
                            IVWrk1,IVWrk2,IVWrk3,IVWrk4);

    DepotCompNr = 0;
    for (i=1; i<=NoOfComponents; i++)
    {
      CompDemand[i] = 0;
      for (k=1; k<=CompsRPtr->LP[i].CFN; k++)
      {
        j = CompsRPtr->LP[i].FAL[k];
        CompNr[j] = i;

        if (j > NoOfCustomers)
        DepotCompNr = i;
        else
        CompDemand[i] += Demand[j];
      }
    }

    for (i=1; i<=NoOfComponents; i++)
    for (j=1; j<=NoOfComponents; j++)
    SMatrix[i][j] = 0.0;

    for (i=1; i<=NoOfCustomers; i++)
    {
      for (k=1; k<=SupportPtr->LP[i].CFN; k++)
      {
        j = SupportPtr->LP[i].FAL[k];
        if (j > i)
        {
          XVal = XMatrix[i][j];

          SMatrix[CompNr[i]][CompNr[j]] += XVal;

          if (CompNr[i] != CompNr[j])
          SMatrix[CompNr[j]][CompNr[i]] += XVal;
        }
      }
    }

    /* Look for a pair of customer components that can be shrunk */

    for (CompI=1; CompI<NoOfComponents; CompI++)
    {
      if (CompI == DepotCompNr) continue;

      for (CompJ=CompI+1; CompJ<=NoOfComponents; CompJ++)
      {
        if (CompJ == DepotCompNr) continue;

        if ((SMatrix[CompI][CompJ] < 0.99) ||
            (SMatrix[CompI][CompJ] > 1.01))
        {
          continue;
        }

        /* We must generate a new 1-edge */

        for (k=1; k<=NoOfComponents; k++)
        {
          if ((k == CompI) || (k == CompJ)) continue;
          /* The depot component is OK to use for k. */

          EdgeSum = SMatrix[CompI][k] + SMatrix[CompJ][k];
          if ((EdgeSum < 0.99) || (EdgeSum > 1.01))
          {
            continue;
          }

          if ((k != DepotCompNr) ||
              ((UseDepotMatch) &&
               ((CompDemand[CompI] + CompDemand[CompJ]) < QMin)))
          {
            ReachAddForwArc(CmprsEdgesRPtr,
                            CompsRPtr->LP[CompI].FAL[1],
                            CompsRPtr->LP[CompJ].FAL[1]);

            ReachAddForwArc(CmprsEdgesRPtr,
                            CompsRPtr->LP[CompJ].FAL[1],
                            CompsRPtr->LP[CompI].FAL[1]);

            ShrinkAgain = 1;
            goto ShrinkAgainLabel;
          }
        } /* k */
      } /* CompJ */
    } /* CompI */

    /* End of pair shrinking */

    /* Look for a triplet of customer components that can be shrunk */

    for (CompI=1; CompI<NoOfComponents-1; CompI++)
    {
      if (CompI == DepotCompNr) continue;

      for (CompJ=CompI+1; CompJ<NoOfComponents; CompJ++)
      {
        if (CompJ == DepotCompNr) continue;
        if (SMatrix[CompI][CompJ] < 0.001) continue;

        for (CompK=CompJ+1; CompK<=NoOfComponents; CompK++)
        {
          if (CompK == DepotCompNr) continue;
          if (SMatrix[CompI][CompK] < 0.001) continue;
          if (SMatrix[CompJ][CompK] < 0.001) continue;

          EdgeSum = SMatrix[CompI][CompJ] +
                    SMatrix[CompI][CompK] +
                    SMatrix[CompJ][CompK];

          if ((EdgeSum < 1.99) || (EdgeSum > 2.01))
          {
            continue;
          }

          /* We must generate a new 1-edge */

          for (k=1; k<=NoOfComponents; k++)
          {
            if ((k == CompI) || (k == CompJ) || (k == CompK)) continue;
            /* The depot component is OK to use for k. */

            EdgeSum = SMatrix[CompI][k] +
                      SMatrix[CompJ][k] +
                      SMatrix[CompK][k];

            if ((EdgeSum < 0.99) || (EdgeSum > 1.01))
            {
              continue;
            }

            if ((k != DepotCompNr) ||
                ((UseDepotMatch) &&
                 ((CompDemand[CompI] +
                   CompDemand[CompJ] +
                   CompDemand[CompK]) < QMin)))
            {
              ReachAddForwArc(CmprsEdgesRPtr,
                              CompsRPtr->LP[CompI].FAL[1],
                              CompsRPtr->LP[CompJ].FAL[1]);

              ReachAddForwArc(CmprsEdgesRPtr,
                              CompsRPtr->LP[CompJ].FAL[1],
                              CompsRPtr->LP[CompI].FAL[1]);

              ReachAddForwArc(CmprsEdgesRPtr,
                              CompsRPtr->LP[CompI].FAL[1],
                              CompsRPtr->LP[CompK].FAL[1]);

              ReachAddForwArc(CmprsEdgesRPtr,
                              CompsRPtr->LP[CompK].FAL[1],
                              CompsRPtr->LP[CompI].FAL[1]);

              ShrinkAgain = 1;
              goto ShrinkAgainLabel;
            }
          } /* k */
        } /* CompK */
      } /* CompJ */
    } /* CompI */

    ShrinkAgainLabel:
    i = 0;

  } while (ShrinkAgain);

  *ResultNoOfComponents = NoOfComponents;

  /* Setup result information */
  Idx = 0;
  for (i=1; i<=NoOfComponents; i++)
  {
    if (i == DepotCompNr) continue;

    Idx++;
    ReachSetForwList(ResultCompsRPtr,
                     CompsRPtr->LP[i].FAL,
                     Idx,
                     CompsRPtr->LP[i].CFN);
  }

  /* Put the depot in the last component */
  Idx++;
  ReachSetForwList(ResultCompsRPtr,
                   CompsRPtr->LP[DepotCompNr].FAL,
                   Idx,
                   CompsRPtr->LP[DepotCompNr].CFN);


  /* Recompute result data for ResultCompsRPtr */

  for (i=1; i<=NoOfComponents; i++)
  {
    CompDemand[i] = 0;
    for (k=1; k<=ResultCompsRPtr->LP[i].CFN; k++)
    {
      j = ResultCompsRPtr->LP[i].FAL[k];
      CompNr[j] = i;

      if (j <= NoOfCustomers)
      CompDemand[i] += Demand[j];
    }
  }

  for (i=1; i<=NoOfComponents; i++)
  for (j=1; j<=NoOfComponents; j++)
  SMatrix[i][j] = 0.0;

  for (i=1; i<=NoOfComponents; i++)
  CompBoundary[i] = 0.0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (k=1; k<=SupportPtr->LP[i].CFN; k++)
    {
      j = SupportPtr->LP[i].FAL[k];
      if (j > i)
      {
        XVal = XMatrix[i][j];

        SMatrix[CompNr[i]][CompNr[j]] += XVal;

        if (CompNr[i] != CompNr[j])
        SMatrix[CompNr[j]][CompNr[i]] += XVal;

        if (CompNr[i] != CompNr[j])
        {
          CompBoundary[CompNr[i]] += XVal;
          CompBoundary[CompNr[j]] += XVal;
        }
      }
    }
  }

  MemFree(CVWrk1);
  MemFree(IVWrk1);
  MemFree(IVWrk2);
  MemFree(IVWrk3);
  MemFree(IVWrk4);

  ReachFreeMem(&CmprsEdgesRPtr);
  ReachFreeMem(&CompsRPtr);
}

void STRCOMB_MainStrengthenedCombs(ReachPtr SupportPtr, /* Original graph */
                                   int NoOfCustomers,
                                   int CAP,
                                   int *Demand, /* Original demand */
                                   int QMin,
                                   double **XMatrix,
                                   int MaxNoOfCuts,
                                   double *MaxViolation,
                                   CnstrMgrPointer CutsCMP)
{
  char UseDepotMatch;
  char OddTeeth,CandidateComb,UseTwoMatching;
  int i,j,k,t,Idx,HIdx;
  int MaxNoOfHandles,NoOfHandles,HandleNr,NoOfTeeth,MaxNoOfTeeth;
  int NoOfComponents,NoOfSuperNodes;
  int SuperNode;
  int NoOfEdges,EdgeNr;
  int HListSize;
  int Head,Tail;
  int ToothRHS,RHSSum,RHS;
  int IntListSize, ExtListSize;
  int NodeListSize;
  int GeneratedCombs;
  int DemandSum;
  int MinVehicles;
  double ToothLHS,LHSSum,LHS;
  double XVal,XSum,HBoundary,ToothSurplus,ToothSurplusSum;

  char *InHandle;
  char *DepotEdgeBound;
  int *EdgeHead, *EdgeTail, *ToothHead, *ToothTail;
  int *HList, *TList;
  int *Index;
  int *IntList, *ExtList;
  int *NodeList;
  int *CompNr;
  int *SuperDemand;
  double *Value, *ToothXSum;
  double *CompBoundary;

  char **InTooth;
  double **SMatrix;

  ReachPtr HandlesRPtr;
  ReachPtr ToothRPtr;
  ReachPtr CompsRPtr;
  ReachPtr SAdjRPtr;
  CnstrMgrPointer TwoMCutsCMP; /* Two-matching cuts. */

  CompNr = MemGetIV(NoOfCustomers+2);
  SuperDemand = MemGetIV(NoOfCustomers+2);
  CompBoundary = MemGetDV(NoOfCustomers+2);

  SMatrix = MemGetDM(NoOfCustomers+2,NoOfCustomers+2);

  ReachInitMem(&CompsRPtr,NoOfCustomers+2);

  DemandSum = 0;
  for (i=1; i<=NoOfCustomers; i++) DemandSum += Demand[i];

  MinVehicles = 0;
  while ((MinVehicles * CAP) < DemandSum) MinVehicles++;

  XSum = 0.0;
  i = NoOfCustomers+1;
  for (k=1; k<=SupportPtr->LP[i].CFN; k++)
  {
    j = SupportPtr->LP[i].FAL[k];
    XVal = XMatrix[i][j];
    XSum += XVal;
  }

  if (XSum <= (2.0 * MinVehicles) + 0.01)
  UseDepotMatch = 1;
  else
  UseDepotMatch = 0;

  STRCOMB_Shrink(SupportPtr,
                 NoOfCustomers,
                 Demand,
                 QMin,
                 UseDepotMatch,
                 XMatrix,
                 CompNr,
                 SuperDemand,
                 CompBoundary,
                 SMatrix,
                 CompsRPtr,
                 &NoOfComponents);

  NoOfSuperNodes = NoOfComponents - 1;

  NodeList = MemGetIV(NoOfCustomers+2);

  ReachInitMem(&SAdjRPtr,NoOfComponents);

  for (i=1; i<=NoOfComponents; i++)
  {
    NodeListSize=0;
    for (j=1; j<=NoOfComponents; j++)
    {
      if (j == i) continue;
      if (SMatrix[i][j] >= 0.01)
      {
        NodeList[++NodeListSize] = j;
      }
    }

    ReachSetForwList(SAdjRPtr,NodeList,i,NodeListSize);
  }

  InHandle = MemGetCV(NoOfSuperNodes+2);
  HList    = MemGetIV(NoOfSuperNodes+1);
  TList    = MemGetIV(NoOfSuperNodes+1);

  MaxNoOfHandles = 2 * NoOfSuperNodes;
  ReachInitMem(&HandlesRPtr,MaxNoOfHandles);

  TwoMCutsCMP = NULL;
  *MaxViolation = 0.0;

  STRCOMB_GenerateHandlesV2(SAdjRPtr,
                            NoOfSuperNodes,
                            SMatrix,
                            HandlesRPtr,
                            MaxNoOfHandles,
                            &NoOfHandles);

  NoOfEdges = 0;
  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=SAdjRPtr->LP[i].CFN; j++)
    {
      k = SAdjRPtr->LP[i].FAL[j];
      if (k > i) NoOfEdges++;
    }
  }

  MaxNoOfTeeth = NoOfEdges;
  InTooth  = MemGetCM(NoOfSuperNodes+2,MaxNoOfTeeth+1);

  ReachInitMem(&ToothRPtr,MaxNoOfTeeth);

  EdgeHead = MemGetIV(NoOfEdges+1);
  EdgeTail = MemGetIV(NoOfEdges+1);
  Index    = MemGetIV(NoOfEdges+1);
  Value    = MemGetDV(NoOfEdges+1);

  EdgeNr = 0;
  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=SAdjRPtr->LP[i].CFN; j++)
    {
      k = SAdjRPtr->LP[i].FAL[j];
      if (k < i) continue;

      if (k > NoOfSuperNodes)
      { /* k is the depot */
        if ((SuperDemand[i] >= QMin) || (UseDepotMatch == 0))
        {
          /* Don't use this edge as a tooth */
          continue; /* next edge */
        }
      }
      else
      { /* k <= NoOfSuperNodes */
        if ((SuperDemand[i] + SuperDemand[k]) > CAP)
        {
          /* Don't use this edge as a tooth */
          continue; /* next edge */
        }
      }

      EdgeNr++;
      EdgeTail[EdgeNr] = i;
      EdgeHead[EdgeNr] = k;
      Index[EdgeNr]    = EdgeNr;
      Value[EdgeNr]    = SMatrix[i][k];
    }
  }

  NoOfEdges = EdgeNr;

  SortIndexDVDec(Index,Value,NoOfEdges);

  ToothHead = MemGetIV(NoOfEdges+1);
  ToothTail = MemGetIV(NoOfEdges+1);
  ToothXSum = MemGetDV(NoOfEdges+1);
  ToothXSum[0] = 0.0;

  UseTwoMatching=0;
  GeneratedCombs=0;

  do
  {
    for (HandleNr=1; HandleNr<=NoOfHandles; HandleNr++)
    {
      for (i=1; i<=NoOfSuperNodes+1; i++) InHandle[i] = 0;

      if (UseTwoMatching)
      {
        HIdx = HandleNr-1;
        HListSize=0;
        for (i=1; i<=TwoMCutsCMP->CPL[HIdx]->IntListSize; i++)
        {
          j = TwoMCutsCMP->CPL[HIdx]->IntList[i];
          InHandle[j] = 1;
          HList[++HListSize] = j;
        }
      }
      else
      {
        HListSize=0;
        for (i=1; i<=HandlesRPtr->LP[HandleNr].CFN; i++)
        {
          j = HandlesRPtr->LP[HandleNr].FAL[i];
          InHandle[j] = 1;
          HList[++HListSize] = j;
        }
      }

      HBoundary = 0.0;
      for (Idx=1; Idx<=HListSize; Idx++)
      {
        i = HList[Idx];
        for (j=1; j<=SAdjRPtr->LP[i].CFN; j++)
        {
          k = SAdjRPtr->LP[i].FAL[j];
          if (InHandle[k] == 0)
          {
            HBoundary += SMatrix[i][k];
          }
        }
      }

      if (UseTwoMatching==0)
      {
        NoOfTeeth=0;
        OddTeeth=0; /* 1 if NoOfTeeth is Odd, 0 otherwise. */
        for (Idx=1; Idx<=NoOfEdges; Idx++)
        {
          EdgeNr = Index[Idx];
          Tail = EdgeTail[EdgeNr];
          Head = EdgeHead[EdgeNr];
          XVal = Value[EdgeNr];

          if ((NoOfTeeth >= 3) && (OddTeeth) && (XVal <= 0.5))
          {
            break;
          }

          if (InHandle[Tail] != InHandle[Head])
          {
            /* Add tooth to list */
            NoOfTeeth++;
            OddTeeth = !OddTeeth;

            ToothXSum[NoOfTeeth] = ToothXSum[NoOfTeeth-1] + XVal;
            ToothTail[NoOfTeeth] = Tail;
            ToothHead[NoOfTeeth] = Head;

            if ((NoOfTeeth > 3) && (OddTeeth))
            {
              if (ToothXSum[NoOfTeeth] <= (ToothXSum[NoOfTeeth-2] + 1.0))
              {
                NoOfTeeth -= 2;
                break;
              }
            }
          }
        }

        if ((NoOfTeeth > 3) && (!OddTeeth))
        {
          NoOfTeeth--;
        }

      }
      else
      {
        NoOfTeeth = TwoMCutsCMP->CPL[HIdx]->ExtListSize / 2;
      }

      if ((NoOfTeeth >= 3) && (HBoundary < (NoOfTeeth + 1.0)))
      {
        /* Candidate handle ok */

        for (j=1; j<=NoOfSuperNodes+1; j++)
        for (t=1; t<=NoOfTeeth; t++)
        InTooth[j][t] = 0;

        ReachClearForwLists(ToothRPtr);

        for (t=1; t<=NoOfTeeth; t++)
        {
          if (UseTwoMatching==0)
          {
            Tail = ToothTail[t];
            Head = ToothHead[t];
          }
          else
          {
            Tail = TwoMCutsCMP->CPL[HIdx]->ExtList[(2*t)-1];
            Head = TwoMCutsCMP->CPL[HIdx]->ExtList[2*t];
          }

          InTooth[Tail][t] = 1;
          InTooth[Head][t] = 1;

          TList[1] = Tail;
          TList[2] = Head;
          ReachSetForwList(ToothRPtr,TList,t,2);
        }

        for (t=1; t<=NoOfTeeth; t++)
        {
          Tail = ToothRPtr->LP[t].FAL[1];
          Head = ToothRPtr->LP[t].FAL[2];
        }

        ToothSurplusSum = 0.0;
        LHSSum = 0.0;
        RHSSum = 0;

        CandidateComb = 1;

        for (t=1; t<=NoOfTeeth; t++)
        {
          /* Try to improve tooth t */

          Tail = ToothRPtr->LP[t].FAL[1];
          Head = ToothRPtr->LP[t].FAL[2];

          STRCOMB_ExpandToothTwoWays(SAdjRPtr,
                                      NoOfSuperNodes,
                                      NoOfTeeth,
                                      t,
                                      SuperDemand,
                                      CAP,
                                      CompBoundary,
                                      InHandle,
                                      InTooth,
                                      SMatrix,
                                      ToothRPtr,
                                      &ToothLHS,
                                      &ToothRHS);

          LHSSum += ToothLHS;
          RHSSum += ToothRHS;

          ToothSurplus = ToothRHS - ToothLHS;
          ToothSurplusSum += ToothSurplus;

          if ((ToothSurplusSum + (NoOfTeeth - t)) <= (HBoundary - 1.0))
          {
            /* Violation is not possible */
            CandidateComb = 0;
            break;
          }
        }

        if (CandidateComb)
        if ((RHSSum % 2) == 0)
        {
          /* RHSSum is even */
          CandidateComb = 0;
        }

        LHSSum += HBoundary;

        if (CandidateComb)
        if (LHSSum < ((1.0 * RHSSum) + 0.99))
        {
          /* Generate comb inequality in cutset form */

          IntList = MemGetIV(NoOfCustomers+1);
          ExtList = MemGetIV(NoOfTeeth*(NoOfCustomers+1));

          if (UseTwoMatching)
          {
            IntListSize=0;
            for (i=1; i<=TwoMCutsCMP->CPL[HIdx]->IntListSize; i++)
            {
              SuperNode = TwoMCutsCMP->CPL[HIdx]->IntList[i];

              for (j=1; j<=CompsRPtr->LP[SuperNode].CFN; j++)
              {
                k = CompsRPtr->LP[SuperNode].FAL[j];
                IntList[++IntListSize] = k;
              }
            }
          }
          else
          {
            IntListSize=0;
            for (i=1; i<=HandlesRPtr->LP[HandleNr].CFN; i++)
            {
              SuperNode = HandlesRPtr->LP[HandleNr].FAL[i];

              for (j=1; j<=CompsRPtr->LP[SuperNode].CFN; j++)
              {
                k = CompsRPtr->LP[SuperNode].FAL[j];
                IntList[++IntListSize] = k;
              }
            }
          }

          ExtListSize = NoOfTeeth;
          for (t=1; t<=NoOfTeeth; t++)
          { /* Generate tooth nr. t */
            ExtList[t] = ExtListSize+1;
            /* Tooth nr. t starts at index ExtList[t]. */

            for (i=1; i<=ToothRPtr->LP[t].CFN; i++)
            {
              SuperNode = ToothRPtr->LP[t].FAL[i];

              for (j=1; j<=CompsRPtr->LP[SuperNode].CFN; j++)
              {
                k = CompsRPtr->LP[SuperNode].FAL[j];
                ExtList[++ExtListSize] = k;
              }
            }
          }

          STRCOMB_ComputeBoundaryLHS(SupportPtr,
                                     NoOfCustomers,
                                     XMatrix,
                                     NoOfTeeth,
                                     IntList,
                                     IntListSize,
                                     ExtList,
                                     ExtListSize,
                                     &LHS);

          STRCOMB_ComputeStrCombRHS(NoOfCustomers,
                                    NoOfTeeth,
                                    Demand,
                                    CAP,
                                    IntList,
                                    IntListSize,
                                    ExtList,
                                    ExtListSize,
                                    &RHS);

          if (LHS < (RHS + 0.99))
          {
            CMGR_AddExtCnstr(CutsCMP,
                             CMGR_CT_STR_COMB,
                             NoOfTeeth,
                             IntListSize,
                             IntList,
                             ExtListSize,
                             ExtList,
                             RHS+1);

            if ((RHS + 1.0 - LHS) > *MaxViolation)
            {
              *MaxViolation = (RHS + 1.0 - LHS);
            }

            GeneratedCombs++;
          }

          MemFree(IntList);
          MemFree(ExtList);
        }
      }

      if (GeneratedCombs >= MaxNoOfCuts)
      {
        UseTwoMatching = 1;
        break;
      }
    } /* HandleNr */

    if (UseTwoMatching)
    {
      UseTwoMatching=0;
      break;
    }

    if ((GeneratedCombs == 0) && (UseTwoMatching==0))
    {
      UseTwoMatching=1;

      CMGR_CreateCMgr(&TwoMCutsCMP,100);

      DepotEdgeBound = MemGetCV(NoOfSuperNodes+1);

      for (i=1; i<=NoOfSuperNodes; i++)
      {
        if ((SuperDemand[i] < QMin) && (UseDepotMatch))
        DepotEdgeBound[i] = 1;
        else
        DepotEdgeBound[i] = 2;
      }

      TWOMATCH_ExactTwoMatchings(SAdjRPtr,
                                 NoOfSuperNodes,
                                 DepotEdgeBound,
                                 SMatrix,
                                 TwoMCutsCMP);

      NoOfHandles = TwoMCutsCMP->Size; /* 0..Size-1 */

      MemFree(DepotEdgeBound);
    }

  } while (UseTwoMatching);

  MemFreeCM(InTooth,NoOfSuperNodes+2);
  MemFreeDM(SMatrix,NoOfCustomers+2);

  MemFree(InHandle);
  MemFree(EdgeHead);
  MemFree(EdgeTail);
  MemFree(ToothHead);
  MemFree(ToothTail);
  MemFree(HList);
  MemFree(TList);
  MemFree(Index);
  MemFree(Value);
  MemFree(ToothXSum);
  MemFree(NodeList);
  MemFree(CompNr);
  MemFree(SuperDemand);
  MemFree(CompBoundary);

  if (TwoMCutsCMP != NULL)
  CMGR_FreeMemCMgr(&TwoMCutsCMP);

  ReachFreeMem(&HandlesRPtr);
  ReachFreeMem(&ToothRPtr);
  ReachFreeMem(&CompsRPtr);
  ReachFreeMem(&SAdjRPtr);
}

