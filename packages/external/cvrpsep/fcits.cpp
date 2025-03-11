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
#include "binpack.h"
#include "cnstrmgr.h"
#include "tolerances.h"

typedef struct
{
  int A;
  int B;
  double X;
} EdgeRec;
typedef EdgeRec *EdgePtr;

typedef struct
{
  int PartitionSize;
  ReachPtr PartitionPtr;
  int EPtrSize;
  int EPtrDim;
  EdgePtr EPtr;
  int ShrunkEdgeNr;
} TreeSearchRec;
typedef TreeSearchRec *TreeSearchPtr;


void FCITS_ComputeLHS(ReachPtr CompsRPtr,
                      int NoOfComps,
                      int NoOfSuperNodes,
                      double **FlowMatrix,
                      ReachPtr FlowRPtr,
                      double *LHS)
{
  int i,j,k;
  int *NodeLabel;
  double *Boundary;

  /* FlowPtr contains only edges between customers, not to the depot */

  NodeLabel = MemGetIV(NoOfSuperNodes+1);
  Boundary = MemGetDV(NoOfSuperNodes+1);

  for (i=1; i<=NoOfSuperNodes; i++) NodeLabel[i] = 0;
  for (i=0; i<=NoOfSuperNodes; i++) Boundary[i] = 0.0;

  for (i=1; i<=NoOfComps; i++)
  {
    for (j=1; j<=CompsRPtr->LP[i].CFN; j++)
    {
      k = CompsRPtr->LP[i].FAL[j];
      NodeLabel[k] = i;
    }
  }

  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=FlowRPtr->LP[i].CFN; j++)
    {
      k = FlowRPtr->LP[i].FAL[j];
      if (k > i)
      {
        if (NodeLabel[k] != NodeLabel[i])
        {
          Boundary[NodeLabel[k]] += FlowMatrix[i][k];
          //if(FlowMatrix[i][k] > 0)
          // printf("Boundary[%d]: %g i:%d k:%d flow:%g\n",
          //        NodeLabel[k],Boundary[NodeLabel[k]],i,k,FlowMatrix[i][k]);
          Boundary[NodeLabel[i]] += FlowMatrix[i][k];
          //if(FlowMatrix[i][k] > 0)
          // printf("Boundary[%d]: %g i:%d k:%d flow:%g\n",
          //        NodeLabel[i],Boundary[NodeLabel[i]],i,k,FlowMatrix[i][k]);
        }
      }
    }
  }

  Boundary[0] = 0.0;
  for (i=1; i<=NoOfSuperNodes; i++)
  {
    Boundary[0] += FlowMatrix[0][i];
    //if(FlowMatrix[0][i] > 0)
    // printf("Boundary[%d]: %g i:%d k:%d flow:%g\n",
    //        0,Boundary[0],0,i,FlowMatrix[0][i]);
    Boundary[NodeLabel[i]] += FlowMatrix[0][i];
    //if(FlowMatrix[0][i] > 0)
    // printf("Boundary[%d]: %g i:%d k:%d flow:%g\n",
    //        NodeLabel[i],Boundary[NodeLabel[i]],0,i,FlowMatrix[0][i]);    
  }

  *LHS = 0.0;
  for (i=0; i<=NoOfComps; i++){
     *LHS += Boundary[i];
     //printf("*LHS:%g Boundary[%d]:%g\n",*LHS,i,Boundary[i]);
  }

  MemFree(NodeLabel);
  MemFree(Boundary);
}

void FCITS_ComputeRHS(ReachPtr PartitionPtr,
                      int PartitionSize,
                      double *SuperNodeDemand,
                      double CAP,
                      double *RHS)
{
  int i,j,k;
  int FullV,LB,UB;
  int *Bin;
  double *ItemSize;

  ItemSize = MemGetDV(PartitionSize+1);
  Bin = MemGetIV(PartitionSize+1);

  FullV = 0;
  for (i=1; i<=PartitionSize; i++)
  {
    ItemSize[i] = 0;
    for (j=1; j<=PartitionPtr->LP[i].CFN; j++)
    {
      k = PartitionPtr->LP[i].FAL[j];
      ItemSize[i] += SuperNodeDemand[k];
      //printf("k:%d i:%d ItemSize:%d SuperNodeDemand:%d\n",
      //   k,i,ItemSize[i],SuperNodeDemand[k]);
    }

    while (ItemSize[i] > CAP)
    {
      FullV++;
      ItemSize[i] -= CAP;
      //printf("FullV=%d ItemSize=%d\n",FullV,ItemSize[i]);
    }
  }

  SortDVDec(ItemSize,PartitionSize);
  BP_ExactBinPacking(CAP,ItemSize,PartitionSize,&LB,&UB,Bin);

  //MVG: rather than start at FullV=1 here, PartitionSize is added in directly,
  //  so, FullV is the number of additional vehicles required for all partitions
  //but, it is not clear to me why that is 4.0 vs 2.0?
  *RHS = (2.0 * (LB + PartitionSize)) + (4.0 * FullV); 
  //printf("RHS:%g LB:%d PartitionSize:%d FullV:%d\n",*RHS,LB,PartitionSize,FullV);

  MemFree(ItemSize);
  MemFree(Bin);
}

void FCITS_ShrinkPartition(int NoOfSuperNodes,
                           ReachPtr FlowPtr,
                           double **FlowMatrix,
                           ReachPtr PartitionPtr,
                           int PartitionSize,
                           int NodeA,
                           int NodeB,
                           ReachPtr NewPartitionPtr,
                           int *NewPartitionSize)
{
  char *InSet, *CompInSet;
  int *CompNr, *NodeList;
  double *XValue;
  int i,j,k;
  int NewNode,NodeListSize;

  InSet = MemGetCV(NoOfSuperNodes+1);
  CompNr = MemGetIV(NoOfSuperNodes+1);
  NodeList = MemGetIV(NoOfSuperNodes+1);

  CompInSet = MemGetCV(PartitionSize+1);
  XValue = MemGetDV(PartitionSize+1);

  NodeListSize=0;

  for (i=1; i<=NoOfSuperNodes; i++) InSet[i] = 0;

  for (i=1; i<=PartitionPtr->LP[NodeA].CFN; i++)
  {
    j = PartitionPtr->LP[NodeA].FAL[i];
    InSet[j]=1;
    NodeList[++NodeListSize] = j;
  }

  for (i=1; i<=PartitionPtr->LP[NodeB].CFN; i++)
  {
    j = PartitionPtr->LP[NodeB].FAL[i];
    InSet[j]=1;
    NodeList[++NodeListSize] = j;
  }

  for (i=1; i<=PartitionSize; i++)
  {
    for (j=1; j<=PartitionPtr->LP[i].CFN; j++)
    {
      k = PartitionPtr->LP[i].FAL[j];
      CompNr[k] = i;
    }
  }

  for (i=1; i<=PartitionSize; i++) XValue[i] = 0.0;

  for (i=1; i<=PartitionSize; i++) CompInSet[i] = 0;
  CompInSet[NodeA] = 1;
  CompInSet[NodeB] = 1;

  do
  {

    for (i=1; i<=PartitionSize; i++) XValue[i] = 0.0;

    for (i=1; i<=NoOfSuperNodes; i++)
    {
      if (InSet[i] == 0) continue;
      for (j=1; j<=FlowPtr->LP[i].CFN; j++)
      {
        k = FlowPtr->LP[i].FAL[j];
        if (InSet[k] == 0)
        XValue[CompNr[k]] += FlowMatrix[i][k];
      }
    }

    NewNode=0;
    for (i=1; i<=PartitionSize; i++)
    if (XValue[i] >= 0.99)
    {
      NewNode = i;
      break;
    }

    if (NewNode > 0)
    {
      /* Induced shrinking */
      for (i=1; i<=PartitionPtr->LP[NewNode].CFN; i++)
      {
        j = PartitionPtr->LP[NewNode].FAL[i];
        InSet[j]=1;
        NodeList[++NodeListSize] = j;
      }

      CompInSet[NewNode] = 1;
    }

  } while (NewNode > 0);


  ReachClearLists(NewPartitionPtr);

  ReachSetForwList(NewPartitionPtr,NodeList,1,NodeListSize);
  *NewPartitionSize = 1;

  for (i=1; i<=PartitionSize; i++)
  {
    if (CompInSet[i] == 0)
    {
      (*NewPartitionSize)++;
      ReachSetForwList(NewPartitionPtr,
                       PartitionPtr->LP[i].FAL,
                       *NewPartitionSize,
                       PartitionPtr->LP[i].CFN);
    }
  }

  MemFree(InSet);
  MemFree(CompNr);
  MemFree(NodeList);

  MemFree(CompInSet);
  MemFree(XValue);
}

void FCITS_CreateEPtrForPartition(int NoOfSuperNodes,
                                  ReachPtr FlowPtr,
                                  double **FlowMatrix,
                                  ReachPtr PartitionPtr,
                                  int PartitionSize,
                                  int Level,
                                  TreeSearchPtr TreePtr)
{
  int i,j,k;
  int EdgeNr,EdgeDim,TotalEdges;
  double XValue;
  int *NodeLabel;
  int *EdgeIndex;
  double *EdgeValue;
  double **EdgeMatrix;
  EdgePtr TmpEdgePtr;

  NodeLabel = MemGetIV(NoOfSuperNodes+1);
  EdgeMatrix = MemGetDM(PartitionSize+1,PartitionSize+1);

  for (i=0; i<=PartitionSize; i++)
  for (j=0; j<=PartitionSize; j++)
  EdgeMatrix[i][j] = 0.0;

  for (i=1; i<=PartitionSize; i++)
  {
    for (j=1; j<=PartitionPtr->LP[i].CFN; j++)
    {
      k = PartitionPtr->LP[i].FAL[j];
      NodeLabel[k] = i;
    }
  }

  EdgeDim = 0;
  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=FlowPtr->LP[i].CFN; j++)
    {
      k = FlowPtr->LP[i].FAL[j];
      if (k > i)
      {
        if (NodeLabel[k] != NodeLabel[i])
        {
          XValue = FlowMatrix[i][k];
          EdgeMatrix[NodeLabel[i]][NodeLabel[k]] += XValue;
          EdgeMatrix[NodeLabel[k]][NodeLabel[i]] += XValue;
          EdgeDim++;
        }
      }
    }
  }

  TmpEdgePtr = (EdgePtr) MemGet(sizeof(EdgeRec)*(EdgeDim+1));
  EdgeIndex = MemGetIV(EdgeDim+1);
  EdgeValue = MemGetDV(EdgeDim+1);

  EdgeNr = 0;
  for (i=1; i<PartitionSize; i++)
  for (j=i+1; j<=PartitionSize; j++)
  {
    if (EdgeMatrix[i][j] >= 0.001)
    {
      EdgeNr++;

      TmpEdgePtr[EdgeNr].A = i;
      TmpEdgePtr[EdgeNr].B = j;
      TmpEdgePtr[EdgeNr].X = EdgeMatrix[i][j];

      EdgeValue[EdgeNr] = TmpEdgePtr[EdgeNr].X;
    }
  }

  TotalEdges = EdgeNr;

  for (i=1; i<=TotalEdges; i++) EdgeIndex[i] = i;
  SortIndexDVDec(EdgeIndex,EdgeValue,TotalEdges);

  if (TreePtr[Level].EPtrDim < TotalEdges)
  {
    TreePtr[Level].EPtr = (EdgePtr) MemReGet(TreePtr[Level].EPtr,
                                             sizeof(EdgeRec)*(TotalEdges+1));
    TreePtr[Level].EPtrDim = TotalEdges;
  }

  TreePtr[Level].EPtrSize = TotalEdges;

  for (i=1; i<=TotalEdges; i++)
  {
    j = EdgeIndex[i];
    TreePtr[Level].EPtr[i].A = TmpEdgePtr[j].A;
    TreePtr[Level].EPtr[i].B = TmpEdgePtr[j].B;
    TreePtr[Level].EPtr[i].X = TmpEdgePtr[j].X;
  }

  MemFree(NodeLabel);
  MemFreeDM(EdgeMatrix,PartitionSize+1);

  MemFree(TmpEdgePtr);
  MemFree(EdgeIndex);
  MemFree(EdgeValue);
}

void FCITS_CheckForDominance(int NoOfSuperNodes,
                             int CurrentLevel,
                             TreeSearchPtr TreePtr,
                             char *Dominated)
{
  int i,j;
  int Level,EdgeNr;
  int NodeA,NodeB,SuperNodeA,SuperNodeB;
  char *InSet;

  InSet = MemGetCV(NoOfSuperNodes+1);
  for (i=1; i<=NoOfSuperNodes; i++) InSet[i] = 0;

  *Dominated = 0;

  /* Check only at index 1 in PartitionPtr. */
  for (i=1; i<=TreePtr[CurrentLevel].PartitionPtr->LP[1].CFN; i++)
  {
    j = TreePtr[CurrentLevel].PartitionPtr->LP[1].FAL[i];
    InSet[j] = 1;
  }

  for (Level=CurrentLevel-1; Level>=0; Level--)
  {
    for (EdgeNr=1; EdgeNr<TreePtr[Level].ShrunkEdgeNr; EdgeNr++)
    {
      NodeA = TreePtr[Level].EPtr[EdgeNr].A;
      NodeB = TreePtr[Level].EPtr[EdgeNr].B;

      SuperNodeA = TreePtr[Level].PartitionPtr->LP[NodeA].FAL[1];
      SuperNodeB = TreePtr[Level].PartitionPtr->LP[NodeB].FAL[1];

      if ((InSet[SuperNodeA]) && (InSet[SuperNodeB]))
      {
        *Dominated = 1;
        goto EndOfDominanceCheck;
      }
    }
  }

  EndOfDominanceCheck:

  MemFree(InSet);
}

void FCITS_TreeSearch(int NoOfSuperNodes,
                      double CAP,
                      double *SuperNodeDemand,
                      double **FlowMatrix,
                      ReachPtr FlowRPtr,
                      int MaxCuts,
                      int MaxFCITSLoops,
                      int *GeneratedCuts,
                      double  EpsViolation,
                      double *MaxViolation,
                      double *CutsRHS,
                      ReachPtr CutsRPtr)
{
  char Backtrack;
  char BacktrackOnSlack;
  char BacktrackOnViolation;
  char Dominated;
  int i,j,k;
  int EdgeNr;
  int MaxLevel,TotalEdges,ListSize,CutNr;
  int Level,TSLoopNr;
  double LHS,RHS,Violation,Eps;
  int *EdgeIndex, *SuperNodeList, *SetSize;
  double *EdgeValue;
  double *LevelViolation;
  EdgePtr TmpEdgePtr;
  TreeSearchPtr TreePtr;

  CutNr = 0;
  Eps = EpsViolation;
  MaxLevel = NoOfSuperNodes-2;

  TotalEdges = 0;
  for (i=1; i<=NoOfSuperNodes; i++) TotalEdges += FlowRPtr->LP[i].CFN;
  TotalEdges /= 2;

  /* Prepare the partition at the root node (level 0) */

  SuperNodeList = MemGetIV(NoOfSuperNodes+1);
  SetSize = MemGetIV(NoOfSuperNodes+1);

  LevelViolation = MemGetDV(MaxLevel+1);

  TmpEdgePtr = (EdgePtr) MemGet(sizeof(EdgeRec)*(TotalEdges+1));
  EdgeIndex = MemGetIV(TotalEdges+1);
  EdgeValue = MemGetDV(TotalEdges+1);

  TreePtr = (TreeSearchPtr) MemGet(sizeof(TreeSearchRec)*(MaxLevel+1));

  TreePtr[0].EPtr = (EdgePtr) MemGet(sizeof(EdgeRec)*(TotalEdges+1));
  TreePtr[0].EPtrSize = TotalEdges;
  TreePtr[0].EPtrDim = TotalEdges;

  for (i=1; i<=MaxLevel; i++)
  {
    TreePtr[i].EPtrSize = 0;
    TreePtr[i].EPtrDim = 0;
    TreePtr[i].EPtr = NULL;
  }

  for (i=0; i<=MaxLevel; i++)
  ReachInitMem(&(TreePtr[i].PartitionPtr),NoOfSuperNodes);

  EdgeNr = 0;
  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=FlowRPtr->LP[i].CFN; j++)
    {
      k = FlowRPtr->LP[i].FAL[j];
      if (k > i)
      {
        EdgeNr++;

        TmpEdgePtr[EdgeNr].A = i;
        TmpEdgePtr[EdgeNr].B = k;
        TmpEdgePtr[EdgeNr].X = FlowMatrix[i][k];

        EdgeValue[EdgeNr] = TmpEdgePtr[EdgeNr].X;
      }
    }
  }

  for (i=1; i<=TotalEdges; i++) EdgeIndex[i] = i;
  SortIndexDVDec(EdgeIndex,EdgeValue,TotalEdges);

  for (i=1; i<=TotalEdges; i++)
  {
    j = EdgeIndex[i];
    TreePtr[0].EPtr[i].A = TmpEdgePtr[j].A;
    TreePtr[0].EPtr[i].B = TmpEdgePtr[j].B;
    TreePtr[0].EPtr[i].X = TmpEdgePtr[j].X;
  }

  for (i=1; i<=NoOfSuperNodes; i++)
  {
    SuperNodeList[1] = i;
    ListSize = 1;
    ReachSetForwList(TreePtr[0].PartitionPtr,SuperNodeList,i,ListSize);
  }

  FCITS_ComputeLHS(TreePtr[0].PartitionPtr,
                   NoOfSuperNodes,
                   NoOfSuperNodes,
                   FlowMatrix,
                   FlowRPtr,
                   &LHS);

  FCITS_ComputeRHS(TreePtr[0].PartitionPtr,
                   NoOfSuperNodes,
                   SuperNodeDemand,
                   CAP,
                   &RHS);

  Violation = RHS - LHS;

  LevelViolation[0] = Violation;

  if (Violation >= Eps)
  {
    ListSize = 0;
    //printf("Violation:%g CutNr:%d\n",Violation,CutNr);

    for (i=1; i<=NoOfSuperNodes; i++)
    {
      SetSize[i] = 1;
      SuperNodeList[++ListSize] = i;
    }

    CutNr++;
    ReachSetForwList(CutsRPtr,SuperNodeList,CutNr,ListSize);
    ReachSetBackwList(CutsRPtr,SetSize,CutNr,NoOfSuperNodes);

    CutsRHS[CutNr] = RHS;

    if (Violation > *MaxViolation)
    *MaxViolation = Violation;

    goto EndOfTreeSearch;
  }

  /* Begin tree search */

  TreePtr[0].PartitionSize = NoOfSuperNodes;
  TreePtr[0].EPtrSize = TotalEdges;
  TreePtr[0].EPtrDim = TotalEdges;
  TreePtr[0].ShrunkEdgeNr = 0;

  Level=0;
  TSLoopNr=0;
  BacktrackOnSlack=0;
  BacktrackOnViolation=0;
  Dominated=0;

  do
  {
    TSLoopNr++;

    do
    {
      /* Find next level from which to shrink and the edge to shrink
         from this level */

      Backtrack=0;

      if (BacktrackOnSlack==1)
      {
        Backtrack=1;
        BacktrackOnSlack=0;
      }
      else
      if (BacktrackOnViolation==1)
      {
        Backtrack=1;
        BacktrackOnViolation=0;
      }
      else
      if (Level >= MaxLevel)
      {
        Backtrack=1;
      }
      else
      if (TreePtr[Level].PartitionSize < 2)
      {
        Backtrack=1;
      }
      else
      if (TreePtr[Level].ShrunkEdgeNr >= TreePtr[Level].EPtrSize)
      {
        Backtrack=1;
      }
      else
      if (Dominated==1)
      {
        Backtrack=1;
        Dominated=0;
      }

      if (Backtrack==1)
      {
        Level--;
        if (Level < 0)
        goto EndOfTreeSearch;
      }
    } while (Backtrack==1);

    if (Level >= MaxLevel)
    goto EndOfTreeSearch;

    (TreePtr[Level].ShrunkEdgeNr)++;

    j = TreePtr[Level].ShrunkEdgeNr;

    FCITS_ShrinkPartition(NoOfSuperNodes,
                          FlowRPtr,
                          FlowMatrix,
                          TreePtr[Level].PartitionPtr,
                          TreePtr[Level].PartitionSize,
                          TreePtr[Level].EPtr[j].A,
                          TreePtr[Level].EPtr[j].B,
                          TreePtr[Level+1].PartitionPtr,
                          &(TreePtr[Level+1].PartitionSize));

    TreePtr[Level+1].ShrunkEdgeNr = 0;

    /* Proceed to next level */
    Level++;

    if (TreePtr[Level].PartitionSize < 2)
    {
      continue;
    }

    /* Check for dominance */

    FCITS_CheckForDominance(NoOfSuperNodes,
                            Level,
                            TreePtr,
                            &Dominated);

    if (Dominated)
    {
      continue;
    }

    FCITS_CreateEPtrForPartition(NoOfSuperNodes,
                                 FlowRPtr,
                                 FlowMatrix,
                                 TreePtr[Level].PartitionPtr,
                                 TreePtr[Level].PartitionSize,
                                 Level,
                                 TreePtr);

    FCITS_ComputeLHS(TreePtr[Level].PartitionPtr,
                     TreePtr[Level].PartitionSize,
                     NoOfSuperNodes,
                     FlowMatrix,
                     FlowRPtr,
                     &LHS);

    FCITS_ComputeRHS(TreePtr[Level].PartitionPtr,
                     TreePtr[Level].PartitionSize,
                     SuperNodeDemand,
                     CAP,
                     &RHS);

    Violation = RHS - LHS;

    LevelViolation[Level] = Violation;

    if (Violation >= Eps)
       {
          //printf("Violation:%g\n",Violation);
      ListSize = 0;
      for (i=1; i<=TreePtr[Level].PartitionSize; i++)
      {
        SetSize[i] = TreePtr[Level].PartitionPtr->LP[i].CFN;
        for (j=1; j<=SetSize[i]; j++)
        {
          k = TreePtr[Level].PartitionPtr->LP[i].FAL[j];
          SuperNodeList[++ListSize] = k;
        }
      }

      CutNr++;
      ReachSetForwList(CutsRPtr,SuperNodeList,CutNr,ListSize);
      ReachSetBackwList(CutsRPtr,SetSize,CutNr,
                        TreePtr[Level].PartitionSize);

      CutsRHS[CutNr] = RHS;

      if (Violation > *MaxViolation)
      *MaxViolation = Violation;

      if (CutNr >= MaxCuts) goto EndOfTreeSearch;

      BacktrackOnViolation = 1;
    }
    else
    if (Violation <= -2.0)
    {
      BacktrackOnSlack=1;
    }

  } while (TSLoopNr < MaxFCITSLoops);

  EndOfTreeSearch:

  *GeneratedCuts = CutNr;

  MemFree(SuperNodeList);
  MemFree(SetSize);

  MemFree(LevelViolation);

  MemFree(TmpEdgePtr);
  MemFree(EdgeIndex);
  MemFree(EdgeValue);

  for (i=0; i<=MaxLevel; i++)
  {
    if (TreePtr[i].EPtr != NULL)
    MemFree(TreePtr[i].EPtr);

    ReachFreeMem(&(TreePtr[i].PartitionPtr));
  }

  MemFree(TreePtr);
}

void FCITS_ComputeFlowMatrix(ReachPtr SupportPtr,
                             int NoOfCustomers,
                             double **XMatrix,
                             ReachPtr SuperNodesRPtr,
                             int NoOfSuperNodes,
                             double **FlowMatrix)
{
  int i,j,k;
  int *NodeLabel;

  NodeLabel = MemGetIV(NoOfCustomers+2);
  for (i=1; i<=NoOfCustomers+1; i++) NodeLabel[i] = 0;

  for (i=1; i<=NoOfSuperNodes; i++)
  {
    for (j=1; j<=SuperNodesRPtr->LP[i].CFN; j++)
    {
      k = SuperNodesRPtr->LP[i].FAL[j];
      NodeLabel[k] = i;
    }
  }

  for (i=0; i<=NoOfSuperNodes; i++)
  for (j=0; j<=NoOfSuperNodes; j++)
  FlowMatrix[i][j] = 0.0;

  for (i=1; i<=NoOfCustomers; i++)
  {
    for (j=1; j<=SupportPtr->LP[i].CFN; j++)
    {
      k = SupportPtr->LP[i].FAL[j];
      if (k > i)
      {
        if (NodeLabel[k] != NodeLabel[i])
        {
          FlowMatrix[NodeLabel[i]][NodeLabel[k]] += XMatrix[i][k];
          FlowMatrix[NodeLabel[k]][NodeLabel[i]] += XMatrix[i][k];

          FlowMatrix[NodeLabel[i]][NodeLabel[i]] += XMatrix[i][k];
          FlowMatrix[NodeLabel[k]][NodeLabel[k]] += XMatrix[i][k];
        }
      }
    }
  }

  MemFree(NodeLabel);
}


void FCITS_MainCutGen(ReachPtr SupportPtr,
                      int NoOfCustomers,
                      const double *Demand,
                      double CAP,
                      double **XMatrix,
                      ReachPtr InitSuperNodesRPtr,
                      ReachPtr InitSAdjRPtr,
                      double *InitSuperDemand,
                      int InitShrunkGraphCustNodes,
                      int MaxFCITSLoops,
                      int MaxNoOfCuts,
                      double  EpsViolation,
                      double *MaxViolation,
                      int *NoOfGeneratedCuts,
                      CnstrMgrPointer CutsCMP)
{
  int i,j,k,Idx,MinIdx,MaxIdx;
  int CompNr;
  double CompDemand;
  int CustNodes,DepotDegree,NoOfComponents;
  int WrkNoOfSuperNodes;
  int NodeListSize;
  int MaxCuts,GeneratedCuts;
  int NoOfSets,SetNr,CutNr;
  int GlobalCutNr;
  double MaxFCIViolation;

  char *CVWrk1;
  int *IVWrk1, *IVWrk2, *IVWrk3, *IVWrk4;

  char *InitSuperNodeInComp;
  double *WrkSuperNodeDemand;
  int *NodeList;
  int *SetSize;
  double *CutsRHS;
  double **FlowMatrix;

  ReachPtr CompsRPtr;
  ReachPtr WrkSuperNodesRPtr;
  ReachPtr FlowRPtr;
  ReachPtr CutsRPtr;

  MaxCuts = MaxNoOfCuts;
  ReachInitMem(&CutsRPtr,MaxCuts);
  CutsRHS = MemGetDV(MaxCuts+1);

  GlobalCutNr = 0;
  MaxFCIViolation = 0.0;

  CustNodes = InitShrunkGraphCustNodes;
  CVWrk1 = MemGetCV(CustNodes+2);
  IVWrk1 = MemGetIV(CustNodes+2);
  IVWrk2 = MemGetIV(CustNodes+2);
  IVWrk3 = MemGetIV(CustNodes+2);
  IVWrk4 = MemGetIV(CustNodes+2);

  InitSuperNodeInComp = MemGetCV(CustNodes+2);
  WrkSuperNodeDemand = MemGetDV(CustNodes+2);

  NodeList = MemGetIV(NoOfCustomers+1);
  SetSize = MemGetIV(CustNodes+2);

  FlowMatrix = MemGetDM(CustNodes+1,CustNodes+1);

  ReachInitMem(&CompsRPtr,CustNodes+1);
  ReachInitMem(&WrkSuperNodesRPtr,CustNodes);
  ReachInitMem(&FlowRPtr,CustNodes);

  DepotDegree = InitSAdjRPtr->LP[CustNodes+1].CFN;
  InitSAdjRPtr->LP[CustNodes+1].CFN = 0;
  ComputeStrongComponents(InitSAdjRPtr,CompsRPtr,
                          &NoOfComponents,CustNodes+1,
                          CVWrk1,
                          IVWrk1,IVWrk2,IVWrk3,IVWrk4);
  InitSAdjRPtr->LP[CustNodes+1].CFN = DepotDegree;

  for (CompNr=1; CompNr<=NoOfComponents; CompNr++)
  {
    if (MaxCuts <= 0) break;
    if (CompsRPtr->LP[CompNr].CFN == 1) continue;

    CompDemand = 0;
    for (i=1; i<=CompsRPtr->LP[CompNr].CFN; i++)
    {
      j = CompsRPtr->LP[CompNr].FAL[i];
      CompDemand += InitSuperDemand[j];
    }

    if (CompDemand <= CAP) continue;

    for (i=1; i<=InitShrunkGraphCustNodes; i++)
    InitSuperNodeInComp[i] = 0;

    for (i=1; i<=CompsRPtr->LP[CompNr].CFN; i++)
    {
      j = CompsRPtr->LP[CompNr].FAL[i];
      InitSuperNodeInComp[j] = 1;
    }

    j=0;
    ReachClearLists(WrkSuperNodesRPtr);
    for (i=1; i<=InitShrunkGraphCustNodes; i++)
    {
      if (InitSuperNodeInComp[i] == 1)
      {
        j++;
        ReachSetForwList(WrkSuperNodesRPtr,
                         InitSuperNodesRPtr->LP[i].FAL,
                         j,
                         InitSuperNodesRPtr->LP[i].CFN);
      }
    }

    WrkNoOfSuperNodes = j;

    FCITS_ComputeFlowMatrix(SupportPtr,
                            NoOfCustomers,
                            XMatrix,
                            WrkSuperNodesRPtr,
                            WrkNoOfSuperNodes,
                            FlowMatrix);

    for (i=1; i<=WrkNoOfSuperNodes; i++)
    {
      NodeListSize = 0;
      for (j=1; j<=WrkNoOfSuperNodes; j++)
      if (j != i)
      if (FlowMatrix[i][j] >= EpsZero)
      {
        NodeList[++NodeListSize] = j;
      }

      if (NodeListSize > 0)
      {
        ReachSetForwList(FlowRPtr,
                         NodeList,
                         i,
                         NodeListSize);
      }
    }

    for (i=1; i<=InitShrunkGraphCustNodes; i++)
    WrkSuperNodeDemand[i] = 0;

    for (k=1; k<=WrkNoOfSuperNodes; k++)
    {
      for (i=1; i<=WrkSuperNodesRPtr->LP[k].CFN; i++)
      {
        j = WrkSuperNodesRPtr->LP[k].FAL[i];
        WrkSuperNodeDemand[k] += Demand[j];
      }
    }

    FCITS_TreeSearch(WrkNoOfSuperNodes,
                     CAP,
                     WrkSuperNodeDemand,
                     FlowMatrix,
                     FlowRPtr,
                     MaxCuts,
                     MaxFCITSLoops,
                     &GeneratedCuts,
                     EpsViolation,
                     &MaxFCIViolation,
                     CutsRHS,
                     CutsRPtr);

    if (GeneratedCuts > 0)
    {
      MaxCuts -= GeneratedCuts;

      for (CutNr=1; CutNr<=GeneratedCuts; CutNr++)
      {
        NodeListSize = 0;
        MaxIdx = 0;
        NoOfSets = CutsRPtr->LP[CutNr].CBN;

        for (SetNr=1; SetNr<=NoOfSets; SetNr++)
        {
          MinIdx = MaxIdx+1;
          MaxIdx = MinIdx + CutsRPtr->LP[CutNr].BAL[SetNr] - 1;

          SetSize[SetNr] = 0;

          for (Idx=MinIdx; Idx<=MaxIdx; Idx++)
          {
            i = CutsRPtr->LP[CutNr].FAL[Idx];
            for (j=1; j<=WrkSuperNodesRPtr->LP[i].CFN; j++)
            {
              k = WrkSuperNodesRPtr->LP[i].FAL[j];
              NodeList[++NodeListSize] = k;
              (SetSize[SetNr])++;
            }
          }
        }

        //printf("Add constraint %d NoOfSets:%d\n",CutNr,NoOfSets);
        CMGR_AddExtCnstr(CutsCMP,
                         CMGR_CT_GENCAP,
                         0,
                         NodeListSize,
                         NodeList,
                         NoOfSets,
                         SetSize,
                         CutsRHS[CutNr]);

        GlobalCutNr++;

      }
    }
  }

  *NoOfGeneratedCuts = GlobalCutNr;
  *MaxViolation = MaxFCIViolation;

  ReachFreeMem(&CutsRPtr);
  MemFree(CutsRHS);

  MemFree(CVWrk1);
  MemFree(IVWrk1);
  MemFree(IVWrk2);
  MemFree(IVWrk3);
  MemFree(IVWrk4);

  MemFree(InitSuperNodeInComp);
  MemFree(WrkSuperNodeDemand);

  MemFree(NodeList);
  MemFree(SetSize);

  MemFreeDM(FlowMatrix,CustNodes+1);

  ReachFreeMem(&CompsRPtr);
  ReachFreeMem(&WrkSuperNodesRPtr);
  ReachFreeMem(&FlowRPtr);
}

