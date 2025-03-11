/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include "memmod.h"
#include "sort.h"
#include "basegrph.h"
#include "compress.h"
#include "cnstrmgr.h"
#include "grsearch.h"
#include "intap.h"
#include "blocks.h"

ReachPtr AllHTourConfigsRPtr = NULL;
int GlobalHTours = 0;
int MaxGlobalHTours = 100;

void NEWHTOUR_CheckIfHandleExists(ReachPtr RPtr,
                                  int RPtrSize, /* #lists in RPtr */
                                  int *HList,
                                  int HListSize,
                                  char *Exists)
{
  int i,j;
  int *SortedList;

  if (RPtrSize < 1)
  {
    *Exists = 0;
    return;
  }

  SortedList = MemGetIV(HListSize+1);
  for (i=1; i<=HListSize; i++) SortedList[i] = HList[i];

  SortIVInc(SortedList,HListSize);

  *Exists = 0;
  for (i=1; i<=RPtrSize; i++)
  {
    if (RPtr->LP[i].CFN != HListSize) continue;

    *Exists = 1;
    for (j=1; j<=HListSize; j++)
    {
      if (RPtr->LP[i].FAL[j] != SortedList[j])
      {
        *Exists = 0;
        break;
      }
    }

    if (*Exists == 1) break;
  }

  MemFree(SortedList);
}

void NEWHTOUR_BuildSets(ReachPtr SupportPtr,
                        int NoOfCustomers,
                        int *Demand, int CAP,
                        double **XMatrix,
                        int NoOfRounds,
                        ReachPtr SetsRPtr,    /* Identified sets */
                        ReachPtr TrySetsRPtr) /* 0=tried before, 1=try this */
{
  char CallBackSets;
  int i,j,k;
  int RoundNr;
  int Source,BestNode,NodeSum;
  int MinCandidateIdx,MaxCandidateIdx;
  int DemandSum;
  int Label;
  int GeneratedSets;
  double BestXVal;
  int *TryVector;
  int *Node, *Pos, *NodeLabel;
  double *XVal;

  TryVector = MemGetIV(NoOfCustomers+1);
  Node = MemGetIV(NoOfCustomers+1);
  Pos = MemGetIV(NoOfCustomers+1);
  NodeLabel = MemGetIV(NoOfCustomers+1);
  XVal = MemGetDV(NoOfCustomers+1);

  for (i=1; i<=NoOfCustomers; i++)
  {
    Node[i] = i;
    Pos[i] = i;
    NodeLabel[i] = 0;
  }

  Label = 0;

  GeneratedSets = 0;

  for (RoundNr=1; RoundNr<=NoOfRounds; RoundNr++)
  for (Source=1; Source<=NoOfCustomers; Source++)
  {
    GRSEARCH_SwapNodesInPos(Node,Pos,1,Pos[Source]);

    if (RoundNr == 1)
    TryVector[1] = 1;
    else
    TryVector[1] = 0;

    DemandSum = Demand[Source];

    MinCandidateIdx = 2;
    MaxCandidateIdx = 1;

    for (j=1; j<=SupportPtr->LP[Source].CFN; j++)
    {
      k = SupportPtr->LP[Source].FAL[j];
      if (k <= NoOfCustomers)
      {
        MaxCandidateIdx++;
        GRSEARCH_SwapNodesInPos(Node,Pos,MaxCandidateIdx,Pos[k]);
        XVal[k] = XMatrix[Source][k];
      }
    }

    NodeSum = Source;
    CallBackSets = 1;
    BestNode = 1;

    while ((MinCandidateIdx <= MaxCandidateIdx) &&
           (BestNode > 0))
    {
      /* The nodes in positions Min...Max are candidates for inclusion. */
      /* Label the nodes that are not feasible for inclusion. */
      Label++;

      /* Call labeling routine. */

      if (CallBackSets)
      {
        GRSEARCH_GetInfeasExt(Pos,
                              MinCandidateIdx,MaxCandidateIdx,
                              NoOfCustomers,
                              NodeSum,
                              SetsRPtr,
                              GeneratedSets,
                              NodeLabel,
                              Label,
                              &CallBackSets);
      }

      BestNode = 0;
      BestXVal = 0.0;

      for (i=MinCandidateIdx; i<=MaxCandidateIdx; i++)
      {
        if (NodeLabel[Node[i]] == Label)
        {
          continue; /* not feasible */
        }

        if ((XVal[Node[i]] > BestXVal) &&
            ((Demand[Node[i]] + DemandSum) <= CAP))
        {
          BestNode = Node[i];
          BestXVal = XVal[Node[i]];
        }
      }

      TryVector[MinCandidateIdx] = 1;

      if (RoundNr > 1)
      if (BestNode == 0)
      { /* Try again, this time allowing copies of earlier sets */
        TryVector[MinCandidateIdx] = 0;

        for (i=MinCandidateIdx; i<=MaxCandidateIdx; i++)
        {
          if ((XVal[Node[i]] > BestXVal) &&
              ((Demand[Node[i]] + DemandSum) <= CAP))
          {
            BestNode = Node[i];
            BestXVal = XVal[Node[i]];
          }
        }
      }

      if (BestNode > 0)
      { /* Include BestNode. */
        GRSEARCH_SwapNodesInPos(Node,Pos,MinCandidateIdx,Pos[BestNode]);
        MinCandidateIdx++;

        NodeSum += BestNode;
        DemandSum += Demand[BestNode];

        /* Update X-values and candidate set. */
        for (j=1; j<=SupportPtr->LP[BestNode].CFN; j++)
        {
          k = SupportPtr->LP[BestNode].FAL[j];
          if (k > NoOfCustomers) continue; /* Depot. */

          if (Pos[k] > MaxCandidateIdx)
          {
            /* k is a new candidate. */
            XVal[k] = XMatrix[BestNode][k];
            MaxCandidateIdx++;
            GRSEARCH_SwapNodesInPos(Node,Pos,MaxCandidateIdx,Pos[k]);
          }
          else
          if (Pos[k] >= MinCandidateIdx)
          {
            /* k is already a candidate */
            XVal[k] += XMatrix[BestNode][k];
          }
          /* Otherwise node k is already in the set */
        }
      }
    }

    GeneratedSets++;
    GRSEARCH_AddSet(SetsRPtr,
                    GeneratedSets,
                    MinCandidateIdx-1,
                    Node,1);
    GRSEARCH_AddSet(TrySetsRPtr,
                    GeneratedSets,
                    MinCandidateIdx-1,
                    TryVector,0);
  }

  MemFree(TryVector);
  MemFree(Node);
  MemFree(Pos);
  MemFree(NodeLabel);
  MemFree(XVal);
}

void NEWHTOUR_QLabelToDepot(ReachPtr SupportPtr,
                            char **AdmissibleEdge,
                            int NoOfCustomers,
                            int *Demand,
                            int TotalDemand,
                            char *CustInSet,
                            int *QLabel,
                            int *NextOnPath)
{
  int i,j,Tmp;
  int LabeledListSize,PermanentLabels;
  int Seed,BestPos,BestLabel,NewLabel;
  char *Labeled;
  int *LabeledList, *Label;

  Labeled     = MemGetCV(NoOfCustomers+1);
  LabeledList = MemGetIV(NoOfCustomers+1);
  Label       = MemGetIV(NoOfCustomers+2);

  LabeledListSize = 0;
  PermanentLabels = 0;

  for (i=1; i<=NoOfCustomers; i++) Labeled[i] = 0;

  Seed = NoOfCustomers+1;
  Label[Seed] = 0;

  do
  {
    /* Label from Seed */
    /* Arc length(i,j) = Demand(j) */

    for (i=1; i<=SupportPtr->LP[Seed].CFN; i++)
    {
      j = SupportPtr->LP[Seed].FAL[i];
      if (AdmissibleEdge[Seed][j] == 0) continue;
      if (j > NoOfCustomers) continue;

      if (CustInSet[Seed] == 0)
      NewLabel = Label[Seed] + Demand[j];
      else
      NewLabel = TotalDemand + Demand[j];

      if ((Labeled[j] == 0) || (NewLabel < Label[j]))
      {
        Label[j] = NewLabel;
        NextOnPath[j] = Seed;

        if (Labeled[j] == 0)
        {
          LabeledList[++LabeledListSize] = j;
          Labeled[j] = 1;
        }
      }
    }

    BestPos=0;
    BestLabel=0;
    for (i=PermanentLabels+1; i<=LabeledListSize; i++)
    {
      j = LabeledList[i];
      if ((BestPos==0) || (Label[j] < BestLabel))
      {
        BestPos = i;
        BestLabel = Label[j];
      }
    }

    if (BestPos > 0)
    {
      Seed = LabeledList[BestPos];

      /* Swap nodes in positions BestPos and PermanentLabels+1 */
      Tmp = LabeledList[BestPos];
      LabeledList[BestPos] = LabeledList[PermanentLabels+1];
      LabeledList[PermanentLabels+1] = Tmp;

      PermanentLabels++;
    }

  } while (BestPos > 0);

  for (i=1; i<=NoOfCustomers; i++)
  QLabel[i] = Label[i] - Demand[i];

  MemFree(Labeled);
  MemFree(LabeledList);
  MemFree(Label);
}

void NEWHTOUR_ComputeBlocks(ReachPtr SupportPtr,
                            char **AdmissibleEdge,
                            int NoOfCustomers,
                            char *NodeInGraph,
                            ReachPtr ResultRPtr,
                            int *NoOfBlocks)
{
  int i,j,k,Idx;
  int NoOfEdges,BlocksCount;
  int ListSize,AdjListSize;
  int *NodeInIndex, *IndexForNode;
  int *AdjList;
  ReachPtr RPtr;
  ReachPtr BlocksRPtr;

  NodeInIndex  = MemGetIV(NoOfCustomers+2);
  IndexForNode = MemGetIV(NoOfCustomers+2);
  AdjList      = MemGetIV(NoOfCustomers+2);

  ListSize=0;
  for (i=1; i<=NoOfCustomers+1; i++)
  {
    if (NodeInGraph[i])
    {
      ListSize++;
      NodeInIndex[ListSize] = i;
      IndexForNode[i] = ListSize;
    }
  }

  ReachInitMem(&RPtr,ListSize);

  NoOfEdges = 0;
  for (Idx=1; Idx<=ListSize; Idx++)
  {
    i = NodeInIndex[Idx];
    AdjListSize = 0;

    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (AdmissibleEdge[i][k] == 0) continue;

      if (NodeInGraph[k])
      {
        AdjList[++AdjListSize] = IndexForNode[k];
      }
    }

    ReachSetForwList(RPtr,AdjList,Idx,AdjListSize);
    NoOfEdges += AdjListSize;
  }

  NoOfEdges /= 2;

  ReachInitMem(&BlocksRPtr,NoOfEdges);

  ComputeBlocks(RPtr,BlocksRPtr,ListSize,&BlocksCount);

  for (i=1; i<=BlocksCount; i++)
  {
    AdjListSize = 0;

    for (j=1; j<=BlocksRPtr->LP[i].CFN; j++)
    {
      k = BlocksRPtr->LP[i].FAL[j];
      AdjList[++AdjListSize] = NodeInIndex[k];
    }

    ReachSetForwList(ResultRPtr,AdjList,i,AdjListSize);
  }

  *NoOfBlocks = BlocksCount;

  MemFree(NodeInIndex);
  MemFree(IndexForNode);
  MemFree(AdjList);

  ReachFreeMem(&RPtr);
  ReachFreeMem(&BlocksRPtr);
}


void NEWHTOUR_CheckBlocks(ReachPtr SupportPtr,
                          char **AdmissibleEdge,
                          int NoOfCustomers,
                          int *Demand,
                          int CAP,
                          double **XMatrix,
                          char *CustInSet,
                          int HandleDemand,
                          int *QToDepot,
                          int Head1,
                          int Head2,
                          ReachPtr BlocksRPtr,
                          int *NextOnPath,
                          double MaxDeleteEdgeSum,
                          int *NoOfEliminatedSets,
                          ReachPtr EliminatedSetsRPtr)
{
  int i,j,k,LoopNr,Origin,Seed,BlockNr,BlockSize;
  int ArtPoint;
  int MaxDemandOnPath;
  int EdgeIdx,EdgesOnPath;
  int NodeListSize;
  int LabeledListSize,PermanentLabels;
  int NewLabel,BestPos,BestLabel;
  int Tmp;
  int NoOfBlocks;
  double XSum;
  char *NodeInGraph;
  char *Labeled;
  int *LabeledList, *Label;
  int *TailOnPath, *HeadOnPath, *BlockOnPath;
  int *NodeList;
  char **InBlock;

  Labeled          = MemGetCV(NoOfCustomers+2);
  LabeledList      = MemGetIV(NoOfCustomers+2);
  Label            = MemGetIV(NoOfCustomers+2);

  NodeInGraph = MemGetCV(NoOfCustomers+2);

  TailOnPath = MemGetIV(NoOfCustomers+2);
  HeadOnPath = MemGetIV(NoOfCustomers+2);
  BlockOnPath = MemGetIV(NoOfCustomers+2);

  NodeList = MemGetIV(NoOfCustomers+2);

  *NoOfEliminatedSets = 0;

  for (LoopNr=1; LoopNr<=2; LoopNr++)
  {
    if (LoopNr==1)
    {
      Origin = Head1;
      MaxDemandOnPath = CAP - (HandleDemand + QToDepot[Head2]);
    }
    else
    {
      Origin = Head2;
      MaxDemandOnPath = CAP - (HandleDemand + QToDepot[Head1]);
    }

    LabeledListSize = 0;
    PermanentLabels = 0;

    for (i=1; i<=NoOfCustomers+1; i++) Labeled[i] = 0;

    Seed = Origin;
    Label[Seed] = 0;

    do
    {
      if (Seed <= NoOfCustomers)
      {
        for (i=1; i<=SupportPtr->LP[Seed].CFN; i++)
        {
          j = SupportPtr->LP[Seed].FAL[i];
          if (AdmissibleEdge[Seed][j] == 0) continue;

          if (j <= NoOfCustomers)
          if (CustInSet[j])
          {
            /* Don't label customers in the handle */
            continue;
          }

          NewLabel = Label[Seed] + Demand[j];

          if ((Labeled[j] == 0) || (NewLabel < Label[j]))
          {
            Label[j] = NewLabel;

            if (Labeled[j] == 0)
            {
              LabeledList[++LabeledListSize] = j;
              Labeled[j] = 1;
            }
          }
        }
      }

      BestPos=0;
      BestLabel=MaxDemandOnPath+1;
      for (i=PermanentLabels+1; i<=LabeledListSize; i++)
      {
        j = LabeledList[i];
        if (Label[j] < BestLabel)
        {
          BestPos = i;
          BestLabel = Label[j];
        }
      }

      if (BestPos > 0)
      {
        Seed = LabeledList[BestPos];

        /* Swap nodes in positions BestPos and PermanentLabels+1 */
        Tmp = LabeledList[BestPos];
        LabeledList[BestPos] = LabeledList[PermanentLabels+1];
        LabeledList[PermanentLabels+1] = Tmp;

        PermanentLabels++;
      }

    } while (BestPos > 0);

    for (i=1; i<=NoOfCustomers+1; i++) NodeInGraph[i] = 0;

    for (i=1; i<=LabeledListSize; i++)
    {
      j = LabeledList[i];
      if (j > NoOfCustomers)
      {
        NodeInGraph[j] = 1; /* Depot */
      }
      else
      if ((Label[j] + QToDepot[j]) <= MaxDemandOnPath)
      {
        NodeInGraph[j] = 1;
      }
    }

    NodeInGraph[Origin] = 1;

    NEWHTOUR_ComputeBlocks(SupportPtr,
                           AdmissibleEdge,
                           NoOfCustomers,
                           NodeInGraph,
                           BlocksRPtr,
                           &NoOfBlocks);

    InBlock = MemGetCM(NoOfBlocks+1,NoOfCustomers+2);
    for (i=1; i<=NoOfBlocks; i++)
    for (j=1; j<=NoOfCustomers+1; j++)
    InBlock[i][j] = 0;

    for (BlockNr=1; BlockNr<=NoOfBlocks; BlockNr++)
    {
      for (j=1; j<=BlocksRPtr->LP[BlockNr].CFN; j++)
      {
        k = BlocksRPtr->LP[BlockNr].FAL[j];
        InBlock[BlockNr][k] = 1;
      }
    }

    EdgesOnPath = 0;
    i = Origin;
    while (i <= NoOfCustomers)
    {
      j = NextOnPath[i];
      EdgesOnPath++;

      TailOnPath[EdgesOnPath] = i;
      HeadOnPath[EdgesOnPath] = j;

      BlockOnPath[EdgesOnPath] = 0;
      for (BlockNr=1; BlockNr<=NoOfBlocks; BlockNr++)
      {
        if ((InBlock[BlockNr][i]) && (InBlock[BlockNr][j]))
        {
          BlockOnPath[EdgesOnPath] = BlockNr;
          break;
        }
      }

      i = j;
    }

    for (EdgeIdx=0; EdgeIdx<=EdgesOnPath; EdgeIdx++)
    {
      /* Check if we change block from edge EdgeIdx to edge EdgeIdx+1 */

      if (EdgeIdx == 0)
      ArtPoint = Origin;
      else
      if (EdgeIdx == EdgesOnPath)
      ArtPoint = NoOfCustomers+1;
      else
      if (BlockOnPath[EdgeIdx] != BlockOnPath[EdgeIdx+1])
      ArtPoint = HeadOnPath[EdgeIdx];
      else
      ArtPoint = 0;

      if (ArtPoint == 0) continue;

      /* Node ArtPoint is an articulation point */

      if (ArtPoint != Origin)
      {
        /* Check edges from entry block to cut point */
        BlockNr = BlockOnPath[EdgeIdx];
        BlockSize = BlocksRPtr->LP[BlockNr].CFN;

        NodeList[1] = ArtPoint;
        NodeListSize = 1;

        XSum = 0.0;

        for (j=1; j<=SupportPtr->LP[ArtPoint].CFN; j++)
        {
          k = SupportPtr->LP[ArtPoint].FAL[j];
          if (AdmissibleEdge[ArtPoint][k] == 0) continue;

          if (InBlock[BlockNr][k])
          {
            NodeList[++NodeListSize] = k;
            XSum += XMatrix[ArtPoint][k];
          }
        }

        if (XSum <= MaxDeleteEdgeSum)
        {
          (*NoOfEliminatedSets)++;

          ReachSetForwList(EliminatedSetsRPtr,
                           NodeList,
                           *NoOfEliminatedSets,
                           NodeListSize);
        }
      }

      if (ArtPoint <= NoOfCustomers)
      {
        /* Check edges from cut point to exit block */
        BlockNr = BlockOnPath[EdgeIdx+1];
        BlockSize = BlocksRPtr->LP[BlockNr].CFN;

        if (BlockSize > 2)
        { /* If BlockSize=2, this edge is checked at the next cut point. */
          /* The other end node is also a cut point. */

          NodeList[1] = ArtPoint;
          NodeListSize = 1;

          XSum = 0.0;

          for (j=1; j<=SupportPtr->LP[ArtPoint].CFN; j++)
          {
            k = SupportPtr->LP[ArtPoint].FAL[j];
            if (AdmissibleEdge[ArtPoint][k] == 0) continue;

            if (InBlock[BlockNr][k])
            {
              NodeList[++NodeListSize] = k;
              XSum += XMatrix[ArtPoint][k];
            }
          }

          if (XSum <= MaxDeleteEdgeSum)
          {
            (*NoOfEliminatedSets)++;

            ReachSetForwList(EliminatedSetsRPtr,
                             NodeList,
                             *NoOfEliminatedSets,
                             NodeListSize);
          }
        }
      }
    }

    MemFreeCM(InBlock,NoOfBlocks+1);
  }

  MemFree(Labeled);
  MemFree(LabeledList);
  MemFree(Label);
  MemFree(NodeInGraph);

  MemFree(TailOnPath);
  MemFree(HeadOnPath);
  MemFree(BlockOnPath);
  MemFree(NodeList);
}

void NEWHTOUR_SolveAP(INTAPPtr AP,
                      ReachPtr SupportPtr,
                      char **AdmissibleEdge,
                      int NoOfCustomers,
                      int *Demand,
                      int CAP,
                      char *CustInSet,
                      int Cust1,
                      int Cust2,
                      int *QOnPaths)
{
  int i,j,k;
  int APValue;
  int InfValue;

  InfValue = CAP+1;

  /* Replace columns <Cust1>,<Cust2> by the depot */

  for (i=1; i<=NoOfCustomers; i++)
  for (j=1; j<=NoOfCustomers; j++)
  AP->c[i][j] = InfValue;

  for (i=1; i<=NoOfCustomers; i++)
  {
    if ((CustInSet[i]) && (i != Cust1) && (i != Cust2)) continue;

    for (k=1; k<=SupportPtr->LP[i].CFN; k++)
    {
      j = SupportPtr->LP[i].FAL[k];
      if (AdmissibleEdge[i][j] == 0) continue;
      if (j > i) /* j may be the depot */
      {
        if (j > NoOfCustomers)
        { /* {i,depot} */
          AP->c[i][Cust1] = 0;
          AP->c[i][Cust2] = 0;
        }
        else
        {
          if ((CustInSet[j]) && (j != Cust1) && (j != Cust2)) continue;

          if ((j != Cust1) && (j != Cust2))
          AP->c[i][j] = Demand[j];

          if ((i != Cust1) && (i != Cust2))
          AP->c[j][i] = Demand[i];
        }
      }
    }
  }

  for (i=1; i<=NoOfCustomers; i++)
  if ((i != Cust1) && (i != Cust2))
  AP->c[i][i] = 0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    if (AP->c[i][i] == 0)
    {
      AP->rall[i] = i;
      AP->call[i] = i;
    }
    else
    {
      AP->rall[i] = 0;
      AP->call[i] = 0;
    }

    AP->u[i] = 0;
    AP->v[i] = 0;
  }

  if (AP->rall[Cust1] == 0) INTAPIterate(AP,Cust1);
  if (AP->rall[Cust2] == 0) INTAPIterate(AP,Cust2);

  APValue = INTAPObjValue(AP);
  *QOnPaths = APValue;
}

void NEWHTOUR_CheckPathsIntersection(int NoOfCustomers,
                                     int Head1,
                                     int Head2,
                                     int *NextOnPath,
                                     char *Intersecting)
{
  int i;
  char *OnPath;

  OnPath = MemGetCV(NoOfCustomers+1);
  for (i=1; i<=NoOfCustomers; i++) OnPath[i] = 0;

  *Intersecting = 0;

  i = Head1;
  while (i <= NoOfCustomers)
  {
    OnPath[i] = 1;
    i = NextOnPath[i];
  }

  i = Head2;
  while (i <= NoOfCustomers)
  {
    if (OnPath[i])
    {
      *Intersecting = 1;
      break;
    }

    i = NextOnPath[i];
  }

  MemFree(OnPath);
}

void NEWHTOUR_ComputeCompletingSet(ReachPtr SupportPtr,
                                   char **AdmissibleEdge,
                                   int NoOfCustomers,
                                   int *Demand,
                                   int TotalDemand,
                                   int CAP,
                                   char *CustInSet,
                                   int HandleDemand,
                                   int Head1,
                                   int Head2,
                                   char **EdgeInF)
{
  int i,j,IdxI,LoopNr,Origin,Seed;
  int MaxDemandOnPath;
  int FirstHead,SecondHead;
  int LabeledListSize,PermanentLabels;
  int NewLabel,BestPos,BestLabel;
  int Tmp;
  char *Labeled;
  int *LabeledList, *Label;
  int *QToDepot, *NextOnPath;

  Labeled     = MemGetCV(NoOfCustomers+2);
  LabeledList = MemGetIV(NoOfCustomers+2);
  Label       = MemGetIV(NoOfCustomers+2);

  for (i=1; i<=NoOfCustomers; i++)
  for (j=i; j<=NoOfCustomers+1; j++)
  {
    EdgeInF[i][j] = 0;
    EdgeInF[j][i] = 0;
  }

  QToDepot   = MemGetIV(NoOfCustomers+1);
  NextOnPath = MemGetIV(NoOfCustomers+1);

  NEWHTOUR_QLabelToDepot(SupportPtr,
                         AdmissibleEdge,
                         NoOfCustomers,
                         Demand,
                         TotalDemand,
                         CustInSet,
                         QToDepot,
                         NextOnPath);

  if (QToDepot[Head1] > QToDepot[Head2])
  {
    FirstHead = Head1;
    SecondHead = Head2;
  }
  else
  {
    FirstHead = Head2;
    SecondHead = Head1;
  }

  for (LoopNr=1; LoopNr<=2; LoopNr++)
  {
    if (LoopNr==1)
    {
      Origin = FirstHead;
      MaxDemandOnPath = CAP - HandleDemand;
    }
    else
    {
      Origin = SecondHead;
      MaxDemandOnPath = CAP - (HandleDemand + QToDepot[FirstHead]);
    }

    if (MaxDemandOnPath < 0) continue;

    LabeledList[0] = Origin;
    LabeledListSize = 0;
    PermanentLabels = 0;

    for (i=1; i<=NoOfCustomers+1; i++) Labeled[i] = 0;

    Seed = Origin;
    Label[Seed] = 0;

    do
    {
      if (Seed <= NoOfCustomers)
      {
        for (i=1; i<=SupportPtr->LP[Seed].CFN; i++)
        {
          j = SupportPtr->LP[Seed].FAL[i];
          if (AdmissibleEdge[Seed][j] == 0) continue;

          if (j <= NoOfCustomers)
          if (CustInSet[j])
          {
            /* Don't label customers in the handle */
            continue;
          }

          NewLabel = Label[Seed] + Demand[j];

          if ((Labeled[j] == 0) || (NewLabel < Label[j]))
          {
            Label[j] = NewLabel;

            if (Labeled[j] == 0)
            {
              LabeledList[++LabeledListSize] = j;
              Labeled[j] = 1;
            }
          }
        }
      }

      BestPos=0;
      BestLabel=MaxDemandOnPath+1;
      for (i=PermanentLabels+1; i<=LabeledListSize; i++)
      {
        j = LabeledList[i];
        if (Label[j] < BestLabel)
        {
          BestPos = i;
          BestLabel = Label[j];
        }
      }

      if (BestPos > 0)
      {
        Seed = LabeledList[BestPos];

        /* Swap nodes in positions BestPos and PermanentLabels+1 */
        Tmp = LabeledList[BestPos];
        LabeledList[BestPos] = LabeledList[PermanentLabels+1];
        LabeledList[PermanentLabels+1] = Tmp;

        PermanentLabels++;
      }

    } while (BestPos > 0);

    for (IdxI=0; IdxI<=PermanentLabels; IdxI++)
    {
      i = LabeledList[IdxI];
      if (i > NoOfCustomers) continue;

      if (AdmissibleEdge[i][NoOfCustomers+1] == 0)
      {
        EdgeInF[i][NoOfCustomers+1] = 1;
        EdgeInF[NoOfCustomers+1][i] = 1;
      }

      for (j=1; j<=NoOfCustomers; j++)
      {
        if (CustInSet[j]) continue;

        if (AdmissibleEdge[i][j] == 0)
        if ((Demand[j] + Label[i] + HandleDemand) <= CAP)
        {
          EdgeInF[i][j] = 1;
          EdgeInF[j][i] = 1;
        }
      }
    }
  }

  MemFree(Labeled);
  MemFree(LabeledList);
  MemFree(Label);

  MemFree(QToDepot);
  MemFree(NextOnPath);
}

void NEWHTOUR_HToursForHandle(ReachPtr SupportPtr,
                              int NoOfCustomers,
                              int *Demand,
                              int CAP,
                              int TotalDemand,
                              double **XMatrix,
                              int *HList,
                              int HListSize,
                              int HDemand,
                              double HBoundary,
                              double *XCustSum,
                              char **AdmissibleEdge,
                              char **EdgeInF,
                              INTAPPtr AP,
                              ReachPtr BlocksRPtr,
                              ReachPtr EliminatedSetsRPtr,
                              int *TailList,
                              int *HeadList,
                              double *CoeffList,
                              int **CutCoeff,
                              CnstrMgrPointer CutsCMP,
                              int MaxGeneratedCuts,
                              char NewHandle,
                              double *MaxViolation)
{
  /* Tries all possible combinations of two edges for this handle */
  char CutFound;
  char ViolationFound,MainViolationFound;
  char ViolationPossible;
  char IntersectingPaths;
  int i,j,k,IdxI,IdxJ;
  int NoOfSpecialPairs,DimSpecialPairs,SpecialPairNr;
  int NoOfEliminatedSets;
  int MinEdgeTail,MinEdgeHead;
  int Head1,Head2;
  int PairQSum;
  int SetNr,NodeListSize;
  int CutListSize;
  int GeneratedCuts;
  double XVal,ViolationUB,DeletedSum,MinEdgeOnPaths;
  double LHS,RHS;
  char *CustInSet;
  int *Mate, *SpecialH1, *SpecialH2, *QToDepot, *NextOnPath, *NodeList;
  double *Alpha, *Score, *MaxDeleteSum;

  CustInSet = MemGetCV(NoOfCustomers+2);
  QToDepot   = MemGetIV(NoOfCustomers+1);
  NextOnPath = MemGetIV(NoOfCustomers+1);

  Mate  = MemGetIV(NoOfCustomers+1);
  Alpha = MemGetDV(NoOfCustomers+1);
  Score = MemGetDV(NoOfCustomers+1);

  DimSpecialPairs = 10;
  SpecialH1 = MemGetIV(DimSpecialPairs+1);
  SpecialH2 = MemGetIV(DimSpecialPairs+1);
  MaxDeleteSum = MemGetDV(DimSpecialPairs+1);

  NodeList = MemGetIV(NoOfCustomers+1);

  /* Let H = {all customers in the set} */
  /* Compute Mate(i) in H, such that x^*(i,Mate(i)) is maximum */

  ViolationPossible = 0;
  GeneratedCuts = 0;

  for (IdxI=1; IdxI<=HListSize; IdxI++)
  {
    i = HList[IdxI];
    Alpha[i] = -1.0;
    Mate[i] = 0;

    for (IdxJ=1; IdxJ<=HListSize; IdxJ++)
    {
      j = HList[IdxJ];
      if (j == i) continue;

      XVal = XMatrix[i][j];
      if (XVal > Alpha[i])
      {
        Alpha[i] = XVal;
        Mate[i] = j;
      }
    }

    Score[i] = 2.0 * (XCustSum[i] - Alpha[i]);

    if ((Score[i] + HBoundary) < 3.98) ViolationPossible = 1;
  }

  if (ViolationPossible == 0) goto EndOfHToursForHandle;

  for (i=1; i<=NoOfCustomers+1; i++) CustInSet[i] = 0;
  for (i=1; i<=HListSize; i++)
  {
    j = HList[i];
    CustInSet[j] = 1;
  }

  NoOfSpecialPairs = 0;
  for (IdxI=1; IdxI<HListSize; IdxI++)
  {
    i = HList[IdxI];

    for (IdxJ=IdxI+1; IdxJ<=HListSize; IdxJ++)
    {
      j = HList[IdxJ];
      if (j == Mate[i]) continue;
      if (Mate[j] == i) continue;

      if ((Mate[j] == Mate[i]) && (HListSize > 3)) continue;
      /* The two edges {i,Mate(i)}, {j,Mate(j)} are node disjoint, */
      /* unless we only have three nodes */

      ViolationUB = 4.0 + (2.0 * XMatrix[i][j]) -
                    HBoundary - Score[i] - Score[j];

      if (ViolationUB >= 0.02)
      {
        NoOfSpecialPairs++;
        if (NoOfSpecialPairs > DimSpecialPairs)
        {
          DimSpecialPairs *= 2;
          SpecialH1 = (int *) MemReGet(SpecialH1,
                                         sizeof(int)*(DimSpecialPairs+1));
          SpecialH2 = (int *) MemReGet(SpecialH2,
                                         sizeof(int)*(DimSpecialPairs+1));
          MaxDeleteSum = (double *) MemReGet(MaxDeleteSum,
                                       sizeof(double)*(DimSpecialPairs+1));
        }

        SpecialH1[NoOfSpecialPairs]  = i;
        SpecialH2[NoOfSpecialPairs]  = j;
        MaxDeleteSum[NoOfSpecialPairs] = (ViolationUB - 0.02) / 2.0;
      }
    }
  }

  /* Consider the configuration for each pair of edges */

  CutFound = 0;

  for (SpecialPairNr=1; SpecialPairNr<=NoOfSpecialPairs; SpecialPairNr++)
  {
    Head1 = SpecialH1[SpecialPairNr];
    Head2 = SpecialH2[SpecialPairNr];

    DeletedSum = 0.0;

    for (i=1; i<=NoOfCustomers; i++)
    {
      for (k=1; k<=SupportPtr->LP[i].CFN; k++)
      {
        j = SupportPtr->LP[i].FAL[k];
        AdmissibleEdge[i][j] = 1;
        AdmissibleEdge[j][i] = 1;
      }
    }

    NEWHTOUR_QLabelToDepot(SupportPtr,
                           AdmissibleEdge,
                           NoOfCustomers,
                           Demand,
                           TotalDemand,
                           CustInSet,
                           QToDepot,
                           NextOnPath);

    PairQSum = QToDepot[Head1] + QToDepot[Head2];

    ViolationFound = 0;
    if ((HDemand + PairQSum) > CAP)
    {
      ViolationFound = 1;
    }
    else
    {
      NEWHTOUR_CheckPathsIntersection(NoOfCustomers,
                                      Head1,
                                      Head2,
                                      NextOnPath,
                                      &IntersectingPaths);

      if (IntersectingPaths)
      {
        NEWHTOUR_SolveAP(AP,
                         SupportPtr,
                         AdmissibleEdge,
                         NoOfCustomers,
                         Demand,
                         CAP,
                         CustInSet,
                         Head1,
                         Head2,
                         &PairQSum);

        if ((HDemand + PairQSum) > CAP)
        {
          ViolationFound = 1;
        }
      }
    }

    MainViolationFound = ViolationFound;
    /* without eliminated edges */

    NoOfEliminatedSets = 0;
    while (ViolationFound == 0)
    {
      /* Compute the minimum x-value among all edges on the two */
      /* paths to the depot. If this exceeds a certain limit, we */
      /* cannot find a violated hypotour inequality by deleting */
      /* an edge. */

      MinEdgeOnPaths = 2.0;

      MinEdgeTail = 0;
      MinEdgeHead = 0;

      if (IntersectingPaths == 0)
      {
        /* Check the two shortest paths */
        i = Head1;
        while (i <= NoOfCustomers)
        {
          j = NextOnPath[i];

          if (XMatrix[i][j] < MinEdgeOnPaths)
          {
            MinEdgeOnPaths = XMatrix[i][j];
            MinEdgeTail = i;
            MinEdgeHead = j;
          }

          i = j;
        }

        i = Head2;
        while (i <= NoOfCustomers)
        {
          j = NextOnPath[i];

          if (XMatrix[i][j] < MinEdgeOnPaths)
          {
            MinEdgeOnPaths = XMatrix[i][j];
            MinEdgeTail = i;
            MinEdgeHead = j;
          }

          i = j;
        }
      }
      else
      {
        /* Use AP paths */

        i = Head1;
        while (i <= NoOfCustomers)
        {
          j = AP->rall[i];
          if ((j == Head1) || (j == Head2)) j = NoOfCustomers+1;

          if (XMatrix[i][j] < MinEdgeOnPaths)
          {
            MinEdgeOnPaths = XMatrix[i][j];
            MinEdgeTail = i;
            MinEdgeHead = j;
          }

          i = j;
        }

        i = Head2;
        while (i <= NoOfCustomers)
        {
          j = AP->rall[i];
          if ((j == Head1) || (j == Head2)) j = NoOfCustomers+1;

          if (XMatrix[i][j] < MinEdgeOnPaths)
          {
            MinEdgeOnPaths = XMatrix[i][j];
            MinEdgeTail = i;
            MinEdgeHead = j;
          }

          i = j;
        }
      }

      if ((MinEdgeOnPaths + DeletedSum) > MaxDeleteSum[SpecialPairNr])
      {
        break;
      }

      ReachClearForwLists(EliminatedSetsRPtr);

      NEWHTOUR_CheckBlocks(SupportPtr,
                           AdmissibleEdge,
                           NoOfCustomers,
                           Demand,
                           CAP,
                           XMatrix,
                           CustInSet,
                           HDemand,
                           QToDepot,
                           Head1,
                           Head2,
                           BlocksRPtr,
                           NextOnPath,
                           MaxDeleteSum[SpecialPairNr] - DeletedSum,
                           &NoOfEliminatedSets,
                           EliminatedSetsRPtr);

      if (NoOfEliminatedSets > 0)
      {
        ViolationFound = 1;
        break;
      }

      if (NoOfEliminatedSets == 0)
      {
        /* Eliminate the x^*-smallest edge on the two paths */

        DeletedSum += MinEdgeOnPaths;

        AdmissibleEdge[MinEdgeTail][MinEdgeHead] = 0;
        AdmissibleEdge[MinEdgeHead][MinEdgeTail] = 0;

        NEWHTOUR_QLabelToDepot(SupportPtr,
                               AdmissibleEdge,
                               NoOfCustomers,
                               Demand,
                               TotalDemand,
                               CustInSet,
                               QToDepot,
                               NextOnPath);

        PairQSum = QToDepot[Head1] + QToDepot[Head2];

        NEWHTOUR_CheckPathsIntersection(NoOfCustomers,
                                        Head1,
                                        Head2,
                                        NextOnPath,
                                        &IntersectingPaths);

        if (IntersectingPaths)
        {
          NEWHTOUR_SolveAP(AP,
                           SupportPtr,
                           AdmissibleEdge,
                           NoOfCustomers,
                           Demand,
                           CAP,
                           CustInSet,
                           Head1,
                           Head2,
                           &PairQSum);
        }

        if ((HDemand + PairQSum) > CAP)
        {
          ViolationFound = 1;
          MainViolationFound = 1;
        }
      }
    } /* while ViolationFound==0 */

    if (ViolationFound)
    {
      /* Generate hypotour inequalities for this handle and edge pair */

      for (SetNr=0; SetNr<=NoOfEliminatedSets; SetNr++)
      {
        /* use SetNr=0, if violation found without eliminating sets */

        if ((SetNr==0) && (MainViolationFound == 0)) continue;

        if (SetNr > 0)
        {
          /* Eliminate additional edges */
          i = EliminatedSetsRPtr->LP[SetNr].FAL[1];

          for (j=2; j<=EliminatedSetsRPtr->LP[SetNr].CFN; j++)
          {
            k = EliminatedSetsRPtr->LP[SetNr].FAL[j];

            AdmissibleEdge[i][k] = 0;
            AdmissibleEdge[k][i] = 0;
          }
        }

        /* Generate inequality */

        NEWHTOUR_ComputeCompletingSet(SupportPtr,
                                      AdmissibleEdge,
                                      NoOfCustomers,
                                      Demand,
                                      TotalDemand,
                                      CAP,
                                      CustInSet,
                                      HDemand, /* incl. Heads */
                                      Head1,
                                      Head2,
                                      EdgeInF);

        CustInSet[Head1] = 0;
        CustInSet[Head2] = 0;

        for (i=1; i<=NoOfCustomers; i++)
        for (j=i; j<=NoOfCustomers+1; j++)
        {
          CutCoeff[i][j] = 0;
          CutCoeff[j][i] = 0;
        }

        LHS = 0.0;
        for (i=1; i<NoOfCustomers; i++)
        for (j=i+1; j<=NoOfCustomers; j++)
        {
          if ((CustInSet[i]) && (CustInSet[j]))
          {
            CutCoeff[i][j] = 1;
            CutCoeff[j][i] = 1;

            LHS += XMatrix[i][j];
          }
        }

        CustInSet[Head1] = 1;
        CustInSet[Head2] = 1;

        CutCoeff[Mate[Head1]][Head1] = 1;
        CutCoeff[Head1][Mate[Head1]] = 1;
        CutCoeff[Mate[Head2]][Head2] = 1;
        CutCoeff[Head2][Mate[Head2]] = 1;

        LHS += XMatrix[Mate[Head1]][Head1];
        LHS += XMatrix[Mate[Head2]][Head2];

        /* Now LHS includes the two special edges */

        for (i=1; i<=NoOfCustomers; i++)
        for (j=i+1; j<=NoOfCustomers+1; j++)
        {
          if (EdgeInF[i][j])
          {
            CutCoeff[i][j] = -1;
            CutCoeff[j][i] = -1;
            LHS -= XMatrix[i][j];
          }
        }

        /* Now LHS includes F */

        CutListSize = 0;
        for (i=1; i<=NoOfCustomers; i++)
        for (j=i+1; j<=NoOfCustomers+1; j++)
        {
          if (CutCoeff[i][j] != 0)
          {
            CutListSize++;

            TailList[CutListSize] = i;
            HeadList[CutListSize] = j;
            CoeffList[CutListSize] = CutCoeff[i][j];
          }
        }

        LHS = 0.0;
        for (i=1; i<=CutListSize; i++)
        LHS += (CoeffList[i] * XMatrix[TailList[i]][HeadList[i]]);

        RHS = (1.0 * HListSize) - 2.0;

        /* Violation = LHS-RHS */
        if ((LHS-RHS) > *MaxViolation) *MaxViolation = LHS-RHS;

        CMGR_AddExplicitCnstr(CutsCMP,
                              CMGR_CT_TWOEDGES_HYPOTOUR,0,
                              CutListSize,
                              TailList,
                              HeadList,
                              CoeffList,
                              RHS);

        CutFound = 1;

        GeneratedCuts++;
        if (GeneratedCuts >= MaxGeneratedCuts)
        {
          goto EndOfSpecialPairs;
        }

        /* Restore eliminated edges */

        if (SetNr > 0)
        {
          i = EliminatedSetsRPtr->LP[SetNr].FAL[1];

          for (j=2; j<=EliminatedSetsRPtr->LP[SetNr].CFN; j++)
          {
            k = EliminatedSetsRPtr->LP[SetNr].FAL[j];

            AdmissibleEdge[i][k] = 1;
            AdmissibleEdge[k][i] = 1;
          }
        }
      }
    }
  }

  EndOfSpecialPairs:

  if ((NewHandle == 1) &&
      (CutFound == 1) &&
      (GlobalHTours < MaxGlobalHTours))
  {
    /* Store the handle for later iterations */
    GlobalHTours++;
    if (GlobalHTours > AllHTourConfigsRPtr->n)
    {
      ReachPtrExpandDim(AllHTourConfigsRPtr,
                        (AllHTourConfigsRPtr->n + NoOfCustomers));
    }

    NodeListSize = 0;
    for (i=1; i<=NoOfCustomers; i++)
    {
      if (CustInSet[i])
      NodeList[++NodeListSize] = i;
    }

    ReachSetForwList(AllHTourConfigsRPtr,
                     NodeList,
                     GlobalHTours,
                     NodeListSize);
  }

  EndOfHToursForHandle:

  MemFree(CustInSet);
  MemFree(QToDepot);
  MemFree(NextOnPath);

  MemFree(Mate);
  MemFree(Alpha);
  MemFree(Score);

  MemFree(SpecialH1);
  MemFree(SpecialH2);
  MemFree(MaxDeleteSum);

  MemFree(NodeList);
}

void NEWHTOUR_ComputeHTours(ReachPtr SupportPtr,
                            int NoOfCustomers,
                            int *Demand,
                            int CAP,
                            double **XMatrix,
                            double **SMatrix,
                            ReachPtr SuperNodesRPtr,
                            ReachPtr SAdjRPtr,
                            int *SuperDemand,
                            int ShrunkGraphCustNodes,
                            CnstrMgrPointer CutsCMP,
                            int MaxHTourCuts,
                            double *MaxHTourViolation)
{
  char HandleExists;
  int i,j,k;
  int Index,SeedNr,SuperSetSize,AddedSuperNode,RemovedSuperNode,DemandSum;
  int CustNr,NoOfRemovedCusts;
  int CustListSize;
  int TotalDemand,TotalNoOfEdges;
  int PreviousHandleNr;
  int CutsBefore,CutsAfter;
  double SBoundary;
  double BestBoundaryFromSeed;
  double MaxViolationForHandle;

  int *CustList;
  int *QToDepot, *NextOnPath;
  int *NodeList;
  double *XCustSum;
  char **EdgeInF;
  char **AdmissibleEdge;
  int **CutCoeff;
  ReachPtr SetsRPtr;
  ReachPtr TrySetsRPtr;

  int CutListSize;
  int *TailList, *HeadList;
  double *CoeffList;

  INTAPPtr AP;

  int NoOfRoundsForSets;

  ReachPtr BlocksRPtr;
  ReachPtr EliminatedSetsRPtr;

  *MaxHTourViolation = 0.0;

  NoOfRoundsForSets = 10;

  if (AllHTourConfigsRPtr == NULL)
  {
    ReachInitMem(&AllHTourConfigsRPtr,NoOfCustomers);
  }

  QToDepot   = MemGetIV(NoOfCustomers+1);
  NextOnPath = MemGetIV(NoOfCustomers+1);

  INTAPInitMem(&AP,NoOfCustomers,NoOfCustomers);

  CustList = MemGetIV(NoOfCustomers+1);
  XCustSum = MemGetDV(NoOfCustomers+1);

  NodeList = MemGetIV(NoOfCustomers+1);

  AdmissibleEdge = MemGetCM(NoOfCustomers+2,NoOfCustomers+2);

  for (i=1; i<=NoOfCustomers; i++)
  for (j=i; j<=NoOfCustomers+1; j++)
  {
    AdmissibleEdge[i][j] = 0;
    AdmissibleEdge[j][i] = 0;
  }

  EdgeInF = MemGetCM(NoOfCustomers+2,NoOfCustomers+2);
  CutCoeff = MemGetIM(NoOfCustomers+2,NoOfCustomers+2);

  CutListSize = ((NoOfCustomers+1) * NoOfCustomers) / 2;
  TailList = MemGetIV(CutListSize+1);
  HeadList = MemGetIV(CutListSize+1);
  CoeffList = MemGetDV(CutListSize+1);

  TotalNoOfEdges = 0;
  for (i=1; i<=NoOfCustomers+1; i++)
  TotalNoOfEdges += SupportPtr->LP[i].CFN;

  TotalNoOfEdges /= 2;

  ReachInitMem(&BlocksRPtr,TotalNoOfEdges);
  ReachInitMem(&SetsRPtr,(ShrunkGraphCustNodes*NoOfRoundsForSets));
  ReachInitMem(&TrySetsRPtr,(ShrunkGraphCustNodes*NoOfRoundsForSets));

  ReachInitMem(&EliminatedSetsRPtr,(NoOfCustomers*4)+1);

  TotalDemand = 0;
  for (i=1; i<=NoOfCustomers; i++)
  TotalDemand += Demand[i];

  /* Try all previously generated handles */

  for (PreviousHandleNr=1;
       PreviousHandleNr<=GlobalHTours;
       PreviousHandleNr++)
  {
    CustListSize = AllHTourConfigsRPtr->LP[PreviousHandleNr].CFN;

    DemandSum = 0;
    for (i=1; i<=CustListSize; i++)
    {
      CustList[i] = AllHTourConfigsRPtr->LP[PreviousHandleNr].FAL[i];
      DemandSum += Demand[CustList[i]];
    }

    for (i=1; i<=NoOfCustomers; i++) XCustSum[i] = 0.0;
    SBoundary = 0.0;

    for (j=1; j<=CustListSize; j++)
    {
      CustNr = CustList[j];

      SBoundary += (2.0 - (2.0 * XCustSum[CustNr]));

      for (k=1; k<=SupportPtr->LP[CustNr].CFN; k++)
      {
        i = SupportPtr->LP[CustNr].FAL[k];
        if (i <= NoOfCustomers)
        {
          XCustSum[i] += XMatrix[CustNr][i];
        }
      }
    }

    if (SBoundary <= 3.98)
    {
      CutsBefore = CutsCMP->Size;
      NEWHTOUR_HToursForHandle(SupportPtr,
                               NoOfCustomers,
                               Demand,
                               CAP,
                               TotalDemand,
                               XMatrix,
                               CustList,
                               CustListSize,
                               DemandSum,
                               SBoundary,
                               XCustSum,
                               AdmissibleEdge,
                               EdgeInF,
                               AP,
                               BlocksRPtr,
                               EliminatedSetsRPtr,
                               TailList,
                               HeadList,
                               CoeffList,
                               CutCoeff,
                               CutsCMP,
                               MaxHTourCuts,
                               0, /* NewHandle == 0 */
                               &MaxViolationForHandle);
      CutsAfter = CutsCMP->Size;
      MaxHTourCuts -= (CutsAfter - CutsBefore);

      if (MaxViolationForHandle > *MaxHTourViolation)
      *MaxHTourViolation = MaxViolationForHandle;

      if (MaxHTourCuts <= 0) goto EndOfAllHandles;
    }
  }

  /* Generate new sets */

  NEWHTOUR_BuildSets(SAdjRPtr,
                     ShrunkGraphCustNodes,
                     SuperDemand,
                     CAP,
                     SMatrix,
                     NoOfRoundsForSets,
                     SetsRPtr,
                     TrySetsRPtr);

  for (Index=1; Index<=(ShrunkGraphCustNodes*NoOfRoundsForSets); Index++)
  {
    SeedNr = Index;

    if (TrySetsRPtr->LP[SeedNr].BAL[1] == 0) continue;
    /* All sets from this seed have been tried before */

    /* Begin with the biggest set from this seed, then iteratively */
    /* remove a supernode */

    DemandSum = 0;
    CustListSize = 0;

    for (i=1; i<=NoOfCustomers; i++) XCustSum[i] = 0.0;
    SBoundary = 0.0;

    BestBoundaryFromSeed = 3.98;

    for (SuperSetSize=1;
         SuperSetSize<=SetsRPtr->LP[SeedNr].CFN;
         SuperSetSize++)
    {
      AddedSuperNode = SetsRPtr->LP[SeedNr].FAL[SuperSetSize];
      DemandSum += SuperDemand[AddedSuperNode];

      /* Expand to customer set */
      for (j=1; j<=SuperNodesRPtr->LP[AddedSuperNode].CFN; j++)
      {
        CustNr = SuperNodesRPtr->LP[AddedSuperNode].FAL[j];
        CustList[++CustListSize] = CustNr;

        SBoundary += (2.0 - (2.0 * XCustSum[CustNr]));

        /* Update X-values */
        for (k=1; k<=SupportPtr->LP[CustNr].CFN; k++)
        {
          i = SupportPtr->LP[CustNr].FAL[k];
          if (i <= NoOfCustomers)
          {
            XCustSum[i] += XMatrix[CustNr][i];
          }
        }
      }
    }

    for (SuperSetSize=SetsRPtr->LP[SeedNr].CFN;
         SuperSetSize>=1;
         SuperSetSize--)
    {
      if (SuperSetSize < SetsRPtr->LP[SeedNr].CFN)
      {
        /* Remove a supernode */
        RemovedSuperNode = SetsRPtr->LP[SeedNr].FAL[SuperSetSize+1];

        DemandSum -= SuperDemand[RemovedSuperNode];
        NoOfRemovedCusts = SuperNodesRPtr->LP[RemovedSuperNode].CFN;

        /* Map to customer set */
        for (j=1; j<=NoOfRemovedCusts; j++)
        {
          CustNr = CustList[CustListSize--];
          SBoundary -= (2.0 - (2.0 * XCustSum[CustNr]));

          /* Update X-values */
          for (k=1; k<=SupportPtr->LP[CustNr].CFN; k++)
          {
            i = SupportPtr->LP[CustNr].FAL[k];
            if (i <= NoOfCustomers)
            {
              XCustSum[i] -= XMatrix[CustNr][i];
            }
          }
        }

      }

      if (TrySetsRPtr->LP[SeedNr].FAL[SuperSetSize] == 0)
      {
        /* Tried before */
        continue;
      }

      if (SBoundary <= (BestBoundaryFromSeed - 0.01))
      {
        BestBoundaryFromSeed = SBoundary;
      }
      else
      {
        /* Dominated by a bigger set from this seed,
            or >= 3.98 <=> max. violation <= 0.01 (in inside form) */
        continue;
      }

      if (CustListSize < 3)
      {
        break;
      }

      NEWHTOUR_CheckIfHandleExists(AllHTourConfigsRPtr,
                                   GlobalHTours,
                                   CustList,
                                   CustListSize,
                                   &HandleExists);

      if (HandleExists == 0)
      {
        CutsBefore = CutsCMP->Size;
        NEWHTOUR_HToursForHandle(SupportPtr,
                                 NoOfCustomers,
                                 Demand,
                                 CAP,
                                 TotalDemand,
                                 XMatrix,
                                 CustList,
                                 CustListSize,
                                 DemandSum,
                                 SBoundary,
                                 XCustSum,
                                 AdmissibleEdge,
                                 EdgeInF,
                                 AP,
                                 BlocksRPtr,
                                 EliminatedSetsRPtr,
                                 TailList,
                                 HeadList,
                                 CoeffList,
                                 CutCoeff,
                                 CutsCMP,
                                 MaxHTourCuts,
                                 1, /* NewHandle */
                                 &MaxViolationForHandle);
        CutsAfter = CutsCMP->Size;
        MaxHTourCuts -= (CutsAfter - CutsBefore);

        if (MaxViolationForHandle > *MaxHTourViolation)
        *MaxHTourViolation = MaxViolationForHandle;

        if (MaxHTourCuts <= 0) goto EndOfAllHandles;
      }
    }
  }

  EndOfAllHandles:

  MemFree(QToDepot);
  MemFree(NextOnPath);

  INTAPFreeMem(&AP);

  MemFree(CustList);
  MemFree(XCustSum);

  MemFree(NodeList);

  MemFreeCM(AdmissibleEdge,NoOfCustomers+2);
  MemFreeCM(EdgeInF,NoOfCustomers+2);
  MemFreeIM(CutCoeff,NoOfCustomers+2);

  MemFree(TailList);
  MemFree(HeadList);
  MemFree(CoeffList);

  ReachFreeMem(&BlocksRPtr);
  ReachFreeMem(&SetsRPtr);
  ReachFreeMem(&TrySetsRPtr);
  ReachFreeMem(&EliminatedSetsRPtr);
}

