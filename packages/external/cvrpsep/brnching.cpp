/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include "memmod.h"
#include "basegrph.h"
#include "sort.h"
#include "cnstrmgr.h"
#include "grsearch.h"
#include "compress.h"
#include "capsep.h"

void BRNCHING_BuildSets(ReachPtr SupportPtr,
                        int NoOfCustomers,
                        int *Demand, int CAP,
                        double **XMatrix,
                        double *NodeBoundary,
                        int *NoOfResultSets,
                        double Target,
                        double *SetBoundary,
                        ReachPtr SetsRPtr,
                        ReachPtr ResultSetsRPtr)
{
  char CallBackSets;
  int i,j,k;
  int Source,BestNode,NodeSum;
  int MinCandidateIdx,MaxCandidateIdx;
  int DemandSum;
  int Label;
  int GeneratedSets;
  double Delta,BestDelta,SBoundary,ThisBoundary;
  int *Node, *Pos, *NodeLabel;
  double *XVal;

  int *NodeList;
  int NodeListSize;

  Node = MemGetIV(NoOfCustomers+1);
  Pos = MemGetIV(NoOfCustomers+1);
  NodeLabel = MemGetIV(NoOfCustomers+1);
  XVal = MemGetDV(NoOfCustomers+1);

  NodeList = MemGetIV(NoOfCustomers+1);

  for (i=1; i<=NoOfCustomers; i++)
  {
    Node[i] = i;
    Pos[i] = i;
    NodeLabel[i] = 0;
  }

  Label = 0;

  GeneratedSets = 0;
  *NoOfResultSets = 0;

  for (Source=1; Source<=NoOfCustomers; Source++)
  {
    GRSEARCH_SwapNodesInPos(Node,Pos,1,Pos[Source]);

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
    SBoundary = NodeBoundary[Source];

    NodeList[1] = Source;
    NodeListSize = 1;

    while ((MinCandidateIdx <= MaxCandidateIdx) &&
           (BestNode > 0))
    {
      Label++;

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
      BestDelta = 2.0;

      for (i=MinCandidateIdx; i<=MaxCandidateIdx; i++)
      {
        if (NodeLabel[Node[i]] == Label)
        {
          continue; /* not feasible */
        }

        ThisBoundary = SBoundary +
                       (NodeBoundary[Node[i]] - (2.0 * XVal[Node[i]]));
        Delta = ThisBoundary - Target;
        if (Delta < 0) Delta = -Delta;

        if ((Delta < BestDelta) &&
            ((Demand[Node[i]] + DemandSum) <= CAP))
        {
          BestNode = Node[i];
          BestDelta = Delta;
        }
      }

      if (BestNode > 0)
      { /* Include BestNode. */
        GRSEARCH_SwapNodesInPos(Node,Pos,MinCandidateIdx,Pos[BestNode]);
        MinCandidateIdx++;

        NodeSum += BestNode;
        DemandSum += Demand[BestNode];

        SBoundary += (NodeBoundary[BestNode] - (2.0 * XVal[BestNode]));

        NodeList[++NodeListSize] = BestNode;
        (*NoOfResultSets)++;
        ReachSetForwList(ResultSetsRPtr,NodeList,*NoOfResultSets,NodeListSize);
        SetBoundary[*NoOfResultSets] = SBoundary;

        for (j=1; j<=SupportPtr->LP[BestNode].CFN; j++)
        {
          k = SupportPtr->LP[BestNode].FAL[j];
          if (k > NoOfCustomers) continue;

          if (Pos[k] > MaxCandidateIdx)
          {
            XVal[k] = XMatrix[BestNode][k];
            MaxCandidateIdx++;
            GRSEARCH_SwapNodesInPos(Node,Pos,MaxCandidateIdx,Pos[k]);
          }
          else
          if (Pos[k] >= MinCandidateIdx)
          {
            XVal[k] += XMatrix[BestNode][k];
          }
        }
      }
    }

    GeneratedSets++;
    GRSEARCH_AddSet(SetsRPtr,
                    GeneratedSets,
                    MinCandidateIdx-1,
                    Node,1);
  }

  MemFree(Node);
  MemFree(Pos);
  MemFree(NodeLabel);
  MemFree(XVal);

  MemFree(NodeList);
}


void BRNCHING_GenerateCandidateSets(ReachPtr SupportPtr,
                                    int NoOfCustomers,
                                    int *Demand,
                                    int CAP,
                                    double BoundaryTarget,
                                    double **XMatrix,
                                    double **SMatrix,
                                    CnstrMgrPointer CMPExistingCuts,
                                    int MaxNoOfCustSets,
                                    ReachPtr *CustSetsRPtr,
                                    int *NoOfCustSets,
                                    double *CustSetBoundary)
{
  int i,j,k,Idx,IdxTail,IdxHead,Tail,Head;
  int MaxSuperNodesInSet,ListDemandSum;
  int NoOfV1Cuts,ShrunkGraphCustNodes;
  int MaxTotalSets,NoOfCandidateSets;
  int TotalEdges;
  int CandidateDemandSum;
  double XVal,XSum;

  int *SuperDemand, *SortedSuperDemand, *CandidateDemand;
  double *SetBoundary, *CandidateBoundary;

  int *CustList;
  int CustListSize,CustNr;

  int CandidateSetNr;

  int *Index;
  double *Value;

  int *CompNr;
  double *CompBoundary;

  ReachPtr SAdjRPtr;
  ReachPtr SuperNodesRPtr;
  ReachPtr V1CutsPtr;

  ReachPtr SetsRPtr;
  ReachPtr ResultSetsRPtr;
  ReachPtr CandidateSetsRPtr;

  /* Generates candidates for cutset branching */

  CompNr = MemGetIV(NoOfCustomers+2);
  CompBoundary = MemGetDV(NoOfCustomers+2);

  TotalEdges = 0;
  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (k > i)
      TotalEdges++;
    }
  }

  /* Generate the shrunk graph for capacity separation */

  CustList = MemGetIV(NoOfCustomers+1);

  V1CutsPtr = NULL;
  CAPSEP_GetOneVehicleCapCuts(CMPExistingCuts,
                              &V1CutsPtr,
                              &NoOfV1Cuts);

  ReachInitMem(&SAdjRPtr,NoOfCustomers+1);
  ReachInitMem(&SuperNodesRPtr,NoOfCustomers+1);

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

  for (i=1; i<=ShrunkGraphCustNodes+1; i++)
  {
    for (k=1; k<=SuperNodesRPtr->LP[i].CFN; k++)
    {
      j = SuperNodesRPtr->LP[i].FAL[k];
      CompNr[j] = i;
    }
  }

  for (i=1; i<=ShrunkGraphCustNodes+1; i++)
  CompBoundary[i] = 0.0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (k=1; k<=SupportPtr->LP[i].CFN; k++)
    {
      j = SupportPtr->LP[i].FAL[k];
      if (j > i)
      {
        XVal = XMatrix[i][j];

        if (CompNr[i] != CompNr[j])
        {
          CompBoundary[CompNr[i]] += XVal;
          CompBoundary[CompNr[j]] += XVal;
        }
      }
    }
  }

  SuperDemand = MemGetIV(ShrunkGraphCustNodes+1);
  SortedSuperDemand = MemGetIV(ShrunkGraphCustNodes+1);

  for (i=1; i<=ShrunkGraphCustNodes; i++)
  {
    SuperDemand[i] = 0;
    for (j=1; j<=SuperNodesRPtr->LP[i].CFN; j++)
    {
      k = SuperNodesRPtr->LP[i].FAL[j];
      SuperDemand[i] += Demand[k];
    }
  }

  for (i=1; i<=ShrunkGraphCustNodes; i++)
  SortedSuperDemand[i] = SuperDemand[i];

  SortIVInc(SortedSuperDemand,ShrunkGraphCustNodes);

  ListDemandSum = SortedSuperDemand[1];
  MaxSuperNodesInSet = 1;

  for (i=2; i<=ShrunkGraphCustNodes; i++)
  {
    if ((SortedSuperDemand[i] + ListDemandSum) > CAP) break;

    MaxSuperNodesInSet++;
    ListDemandSum += SortedSuperDemand[i];
  }

  MaxTotalSets = (MaxSuperNodesInSet * ShrunkGraphCustNodes) + TotalEdges;

  SetBoundary = MemGetDV(MaxTotalSets);
  CandidateBoundary = MemGetDV(MaxTotalSets);

  CandidateDemand = MemGetIV(MaxTotalSets);

  ReachInitMem(&SetsRPtr,ShrunkGraphCustNodes);
  ReachInitMem(&ResultSetsRPtr,MaxTotalSets);
  ReachInitMem(&CandidateSetsRPtr,MaxTotalSets);

  BRNCHING_BuildSets(SAdjRPtr,
                     ShrunkGraphCustNodes,
                     SuperDemand,
                     CAP,
                     SMatrix,
                     CompBoundary,
                     &NoOfCandidateSets,
                     BoundaryTarget,
                     SetBoundary,
                     SetsRPtr,
                     ResultSetsRPtr);

  /* Expand from supernodes to customer set */
  CandidateSetNr = 0;

  for (i=1; i<=NoOfCandidateSets; i++)
  {
    CustListSize = 0;
    CandidateDemandSum = 0;

    for (j=1; j<=ResultSetsRPtr->LP[i].CFN; j++)
    {
      k = ResultSetsRPtr->LP[i].FAL[j]; /* Supernode nr. k */
      /* Expand supernode k */
      for (Idx=1; Idx<=SuperNodesRPtr->LP[k].CFN; Idx++)
      {
        CustNr = SuperNodesRPtr->LP[k].FAL[Idx];
        CustList[++CustListSize] = CustNr;
        CandidateDemandSum += Demand[CustNr];
      }
    }

    if (CustListSize > 2)
    if ((SetBoundary[i] >= 2.002) && (SetBoundary[i] <= 3.998))
    {
      XSum = 0.0;
      for (IdxTail=1; IdxTail<CustListSize; IdxTail++)
      {
        Tail = CustList[IdxTail];
        for (IdxHead=IdxTail+1; IdxHead<=CustListSize; IdxHead++)
        {
          Head = CustList[IdxHead];
          XSum += XMatrix[Tail][Head];
        }
      }

      if ((XSum <= ((1.0 * CustListSize) - 1.001)) &&
          (XSum >= ((1.0 * CustListSize) - 1.999)))
      {
        /* Add to candidate list */
        CandidateSetNr++;
        ReachSetForwList(CandidateSetsRPtr,
                          CustList,
                          CandidateSetNr,
                          CustListSize);

        CandidateBoundary[CandidateSetNr] = SetBoundary[i];
        CandidateDemand[CandidateSetNr] = CandidateDemandSum;
      }
    }
  }

  /* Add edges between customers in original support graph */
  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if ((k > i) && (k <= NoOfCustomers))
      if ((XMatrix[i][k] >= 0.0001) && (XMatrix[i][k] <= 0.9999))
      {
        CustList[1] = i;
        CustList[2] = k;
        CustListSize = 2;

        CandidateSetNr++;
        ReachSetForwList(CandidateSetsRPtr,
                         CustList,
                         CandidateSetNr,
                         CustListSize);

        CandidateBoundary[CandidateSetNr] = 2.0 * (2.0 - XMatrix[i][k]);
        CandidateDemand[CandidateSetNr] = Demand[i] + Demand[k];
      }
    }
  }

  NoOfCandidateSets = CandidateSetNr;

  Index = MemGetIV(NoOfCandidateSets+1);
  Value = MemGetDV(NoOfCandidateSets+1);

  for (i=1; i<=NoOfCandidateSets; i++)
  {
    Index[i] = i;
    Value[i] = (CandidateBoundary[i] - BoundaryTarget);
    if (Value[i] < 0.0) Value[i] = -Value[i];
    Value[i] += (0.0001); /* For tie-breaking if Value = 0.0 */
    Value[i] = (Value[i] / (1.0 * CandidateDemand[i]));
  }

  SortIndexDVInc(Index,Value,NoOfCandidateSets);

  ReachInitMem(CustSetsRPtr,NoOfCandidateSets);

  CandidateSetNr = 0;
  for (i=1; i<=NoOfCandidateSets; i++)
  {
    j = Index[i];

    CustListSize = 0;
    for (k=1; k<=CandidateSetsRPtr->LP[j].CFN; k++)
    {
      CustNr = CandidateSetsRPtr->LP[j].FAL[k];
      CustList[++CustListSize] = CustNr;
    }

    CandidateSetNr++;
    ReachSetForwList(*CustSetsRPtr,CustList,CandidateSetNr,CustListSize);

    CustSetBoundary[CandidateSetNr] = CandidateBoundary[j];

    if (CandidateSetNr >= MaxNoOfCustSets)
    break;
  }

  *NoOfCustSets = CandidateSetNr;

  MemFree(SuperDemand);
  MemFree(SortedSuperDemand);
  MemFree(SetBoundary);
  MemFree(CandidateBoundary);
  MemFree(CustList);
  MemFree(Index);
  MemFree(Value);

  MemFree(CandidateDemand);

  MemFree(CompNr);
  MemFree(CompBoundary);

  ReachFreeMem(&SAdjRPtr);
  ReachFreeMem(&SuperNodesRPtr);
  ReachFreeMem(&SetsRPtr);
  ReachFreeMem(&ResultSetsRPtr);
  ReachFreeMem(&CandidateSetsRPtr);
}

void BRNCHING_GetCandidateSets(int NoOfCustomers,
                               int *Demand,
                               int CAP,
                               int NoOfEdges,
                               int *EdgeTail,
                               int *EdgeHead,
                               double *EdgeX,
                               CnstrMgrPointer CMPExistingCuts,
                               double BoundaryTarget,
                               int MaxNoOfSets,
                               CnstrMgrPointer SetsCMP)
{
  int i,j,SetNr,NoOfSets;
  double *SetBoundary;
  double **XMatrix;
  double **SMatrix;
  ReachPtr SupportPtr;
  ReachPtr CustSetsRPtr;

  ReachInitMem(&SupportPtr,NoOfCustomers+1);

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

  SetBoundary = MemGetDV(MaxNoOfSets+1);

  BRNCHING_GenerateCandidateSets(SupportPtr,
                                 NoOfCustomers,
                                 Demand,
                                 CAP,
                                 BoundaryTarget,
                                 XMatrix,
                                 SMatrix,
                                 CMPExistingCuts,
                                 MaxNoOfSets,
                                 &CustSetsRPtr,
                                 &NoOfSets,
                                 SetBoundary);

  SetsCMP->Size = 0;
  for (SetNr=1; SetNr<=NoOfSets; SetNr++)
  {
    CMGR_AddCnstr(SetsCMP,
                  0,0,
                  CustSetsRPtr->LP[SetNr].CFN,
                  CustSetsRPtr->LP[SetNr].FAL,
                  SetBoundary[SetNr]);
  }

  MemFree(SetBoundary);

  MemFreeDM(SMatrix,NoOfCustomers+2);
  MemFreeDM(XMatrix,NoOfCustomers+2);

  ReachFreeMem(&SupportPtr);
  ReachFreeMem(&CustSetsRPtr);
}

