/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#include <stdlib.h>
#include <stdio.h>
#include "memmod.h"
#include "basegrph.h"
#include "sort.h"

void BP_FirstFit(double CAP,
                 double *ItemSize,
                 int n,
                 int *Bin,
                 int *NoOfBins)
{
  /* This is the First-Fit heuristic algorithm for the
     Bin Packing Problem. If the items are ordered in nonincreasing order,
     the algorithm is the First-Fit Decreasing algorithm.
     Return values: Bin[i] is the number of the bin into which
     item i is packed. */

  int i,Item,BinNr,UsedBins;
  double *Slack;

  /* printf("BEGIN BP_FirstFit\n"); */

  Slack = MemGetDV(n+1);

  UsedBins = 0;

  for (Item=1; Item<=n; Item++)
  {
    BinNr = 0;
    for (i=1; i<=UsedBins; i++)
    {
      if (Slack[i] >= ItemSize[Item])
      {
        BinNr = i;
        break;
      }
    }

    /* printf("Item %d: BinNr = %d\n",Item,BinNr); */

    if (BinNr == 0)
    {
      UsedBins++;
      Slack[UsedBins] = CAP;
      BinNr = UsedBins;
      /* printf("New bin %d\n",BinNr); */
    }

    Bin[Item] = BinNr;
    Slack[BinNr] -= ItemSize[Item];

    /* printf("Slack(%d) = %d\n",BinNr,Slack[BinNr]); */
  }

  /*
  printf("UsedBins = %d\n",UsedBins);
  printf("Bin =");
  for (i=1; i<=n; i++) printf(" %d",Bin[i]);
  printf("\n");
  printf("Slack =");
  for (i=1; i<=UsedBins; i++) printf(" %d",Slack[i]);
  printf("\n");
  */

  *NoOfBins = UsedBins;

  MemFree(Slack);

  /* printf("END BP_FirstFit(Bins=%d)\n",*NoOfBins); */
}

void BP_ModifiedFirstFit(double CAP,
                         double *ItemSize,
                         int n,
                         int *Bin,
                         int *NoOfBins)
{
  /* This is the First-Fit heuristic algorithm for the
     Bin Packing Problem. If the items are ordered in nonincreasing order,
     the algorithm is the First-Fit Decreasing algorithm.
     Return values: Bin[i] is the number of the bin into which
     item i is packed. If Bin[i] > 0 on input, then that allocation is
     fixed. */

  int i,Item,BinNr,UsedBins;
  double *Slack;

  /*
  printf("BEGIN BP_ModifiedFirstFit: CAP=%d, n=%d\n",CAP,n);
  printf("ItemSize =");
  for (i=1; i<=n; i++) printf(" %d",ItemSize[i]);
  printf("\n");
  */

  Slack = MemGetDV(n+1);

  UsedBins = 0;

  for (Item=1; Item<=n; Item++)
  if (Bin[Item] > 0)
  {
    if (Bin[Item] > UsedBins)
    {
      for (i=UsedBins+1; i<=Bin[Item]; i++) Slack[i] = CAP;
      UsedBins = Bin[Item];
    }

    Slack[Bin[Item]] -= ItemSize[Item];
  }

  for (Item=1; Item<=n; Item++)
  {
    if (Bin[Item] > 0) continue;

    BinNr = 0;
    for (i=1; i<=UsedBins; i++)
    {
      if (Slack[i] >= ItemSize[Item])
      {
        BinNr = i;
        break;
      }
    }

    /* printf("Item %d: BinNr = %d\n",Item,BinNr); */

    if (BinNr == 0)
    {
      UsedBins++;
      Slack[UsedBins] = CAP;
      BinNr = UsedBins;
      /* printf("New bin %d\n",BinNr); */
    }

    Bin[Item] = BinNr;
    Slack[BinNr] -= ItemSize[Item];

    /* printf("Slack(%d) = %d\n",BinNr,Slack[BinNr]); */
  }

  /*
  printf("UsedBins = %d\n",UsedBins);
  printf("Bin =");
  for (i=1; i<=n; i++) printf(" %d",Bin[i]);
  printf("\n");
  printf("Slack =");
  for (i=1; i<=UsedBins; i++) printf(" %d",Slack[i]);
  printf("\n");
  */

  *NoOfBins = UsedBins;

  MemFree(Slack);

  /*
  printf("FFD: Bin =");
  for (i=1; i<=n; i++) printf(" %d",Bin[i]);
  printf("\n");
  */

  /* printf("END BP_ModifiedFirstFit(Bins=%d)\n",*NoOfBins); */
}



void BP_DominancePacking(double CAP,
                         double *ItemSize,
                         int n,
                         int *Bin,
                         int *NoOfBins)
{
  /* This is the reduction algorithm in Martello & Toth
     ("Knapsack Problems", 1990, pp. 233-234), based on dominance
     between different packings.
       On return, Bin[i] is the bin number in which item i has been packed.
     If Bin[i] = 0, then item i has not been packed into a bin. The
     remaining items are therefore those with Bin[i] = 0.
       NoOfBins is the number of bins that have been optimally packed.

     It is assumed that on input the items are sorted in nonincreasing
     order of size, i.e., ItemSize[i] >= ItemSize[i+1], and that
     all items are <= the capacity CAP.
  */

  int i,j,k,r,s,PrevJ,JStar;
  int UsedBins;
  int BestR,BestS,PrevR,PrevS;
  int ItemsBetweenRS;
  double RemainingWeight,SumPairBetweenRS,BestSum,CapForS,SumSize,LimitSize;

  /*
  printf("BP_DominancePacking: Input:\n");
  printf("CAP=%d, n=%d\n",CAP,n);
  printf("ItemSize =");
  for (i=1; i<=n; i++) printf(" %d",ItemSize[i]);
  printf("\n");
  */

  RemainingWeight = 0;
  UsedBins = 0;
  PrevJ = 0;

  for (i=1; i<=n; i++)
  {
    Bin[i] = 0;
    RemainingWeight += ItemSize[i];
  }

  while ((PrevJ < n) && (RemainingWeight > CAP))
  {
    j = 0;
    for (i=PrevJ+1; i<=n; i++)
    {
      if (Bin[i]==0)
      {
        j = i;
        break;
      }
    }

    if (j == 0) break;
    PrevJ = j;

    LimitSize = CAP - ItemSize[j];

    JStar = 0;
    for (i=1; i<=n; i++)
    if ((i!=j) && (Bin[i]==0))
    {
      if (ItemSize[i] <= LimitSize)
      {
        JStar = i;
        break;
      }
    }

    if (JStar == 0)
    { /* Packing of item j as the only item into a bin */
      UsedBins++;
      Bin[j] = UsedBins;
      RemainingWeight -= ItemSize[j];
    }
    else
    {
      /* Find the maximum number of items (k) that can be packed together
         with item j. */

      k=0;
      SumSize=0;

      if (ItemSize[JStar] != LimitSize)
      {
        for (i=n; i>=JStar; i--)
        if (i!=j)
        {
          if (Bin[i]==0)
          {
            SumSize += ItemSize[i];
            k++;
            if (SumSize > LimitSize)
            {
              k--;
              break;
            }
            if (k == 3) break;
          }
        }
      }

      if ((k==1) || (ItemSize[JStar] == LimitSize))
      {
        /* Packing of items j and JStar as the only two items into a bin */
        UsedBins++;
        Bin[j] = UsedBins;
        RemainingWeight -= ItemSize[j];
        Bin[JStar] = UsedBins;
        RemainingWeight -= ItemSize[JStar];
      }
      else
      if (k==2)
      {
        /* Left pointer is r, right pointer is s. */

        for (s=n; s>=JStar; s--)
        {
          if ((s!=j) && (Bin[s]==0))
          break;
        }

        /* s is the index of the smallest item */

        for (r=JStar; r<s; r++)
        {
          if ((r!=j) &&
              (Bin[r]==0) && ((ItemSize[r] + ItemSize[s]) <= LimitSize))
          break;
        }

        /* r is the index of the biggest item that can be packed
           together with items j and s */

        BestSum = ItemSize[r] + ItemSize[s];
        BestR = r;
        BestS = s;

        PrevR = r;
        PrevS = s;

        /* The previous r-value is included in the interval for r in
           the following loop, as the biggest item s that can be packed
           together with j and r has not (necessarily) been found yet */

        for (r=PrevR; r<PrevS; r++)
        {
          if ((r!=j) && (Bin[r]==0))
          {
            /* Find the biggest item s that can be packed together with
               j and r */

            CapForS = LimitSize - ItemSize[r];

            for (s=PrevS-1; s>r; s--)
            {
              if ((s!=j) && (Bin[s]==0))
              {
                if (ItemSize[s] > CapForS) break;

                if (ItemSize[s] + ItemSize[r] > BestSum)
                {
                  BestSum = ItemSize[s] + ItemSize[r];
                  BestR   = r;
                  BestS   = s;
                }
              }
            }

            PrevS = BestS;
          }
        }

        /* End of search for k=2 */

        if (ItemSize[JStar] >= BestSum)
        {
          /* JStar is better than r+s */
          /* Packing of items j and JStar as the only two items into a bin */

          UsedBins++;
          Bin[j] = UsedBins;
          RemainingWeight -= ItemSize[j];
          Bin[JStar] = UsedBins;
          RemainingWeight -= ItemSize[JStar];
        }
        else
        if (BestR == JStar)
        {
          ItemsBetweenRS = 0;
          SumPairBetweenRS = 0;

          for (i=BestS-1; i>BestR; i--)
          {
            if ((i!=j) && (Bin[i]==0))
            {
              if ((ItemSize[i] != ItemSize[BestR]) &&
                  (ItemSize[i] != ItemSize[BestS]))
              {
                ItemsBetweenRS++;
                SumPairBetweenRS += ItemSize[i];

                if (ItemsBetweenRS == 2) break;
              }
            }
          }

          if ((ItemsBetweenRS <= 1) || (SumPairBetweenRS > LimitSize))
          {
            /* Packing of items j, BestR and BestS as the
               only three items into a bin */

            UsedBins++;
            Bin[j] = UsedBins;
            RemainingWeight -= ItemSize[j];
            Bin[BestR] = UsedBins;
            RemainingWeight -= ItemSize[BestR];
            Bin[BestS] = UsedBins;
            RemainingWeight -= ItemSize[BestS];
          }
        }
      } /* k = 2 */
    }
  }

  if (RemainingWeight <= CAP)
  {
    /* Pack the remaining items into one bin */
    UsedBins++;
    for (i=1; i<=n; i++)
    {
      if (Bin[i] == 0)
      Bin[i] = UsedBins;
    }
  }

  *NoOfBins = UsedBins;

  /*
  printf("NoOfBins = %d\n",*NoOfBins);
  if ((*NoOfBins) > 0)
  {
    printf("Bin =");
    for (i=1; i<=n; i++) printf(" %d",Bin[i]);
    printf("\n");
  }
  */

  /* printf("END BP_DominancePacking: NoOfBins = %d\n",*NoOfBins); */
}


void BP_MTL2(double CAP,
             double *ItemSize,
             int n,
             int *LB)
{
  /* This is the bounding procedure L2 in Martello & Toth
     ("Knapsack Problems", 1990, pp. 231-232).
       It is assumed that on input the items are sorted in nonincreasing
     order of size, i.e., ItemSize[i] >= ItemSize[i+1], and that
     all items are <= the capacity CAP.
  */

   int JStar = -1;
  int i,CJ12,CJ2,JPrime,JDPrime;
  int /*CAPHalf,*/ AddBins;
  int L2;
  double SJ2,SJ3,Limit,AddCAP,AddSizeSum,SJStar,CAPSum,SumAllItems;


  /*
  printf("BP_MTL2: Input:\n");
  printf("CAP=%d, n=%d\n",CAP,n);
  printf("ItemSize =");
  for (i=1; i<=n; i++) printf(" %d",ItemSize[i]);
  printf("\n");
  */

  /*
  CAPHalf = CAP / 2;
  if (ItemSize[n] > CAPHalf)
  {
    *LB = n;
    return;
  }
  */

  if ((ItemSize[n-1] + ItemSize[n]) > CAP)
  {
    *LB = n;
    return;
  }


  for (i=1; i<=n; i++)
  {
    /*
    if (ItemSize[i] < CAPHalf)
    {
      JStar = i;
      break;
    }
    */

    if ((ItemSize[i] * 2) <= CAP)
    {
      JStar = i;
      break;
    }

  }

  /* printf("JStar = %d\n",JStar); */

  if (JStar == 1)
  {
    /* The bound reduces to the trivial lower bound L1 */
    SumAllItems = 0;
    for (i=1; i<=n; i++) SumAllItems += ItemSize[i];

    (*LB) = 1;
    CAPSum = CAP;
    while (CAPSum < SumAllItems)
    {
      (*LB)++;
      CAPSum += CAP;
    }

    /* printf("Trivial bound L2:= %d\n",*LB); */
    return;
  }

  CJ12 = JStar - 1;

  SJStar = 0;
  for (i=JStar; i<=n; i++) SJStar += ItemSize[i];

  JPrime = 0;
  Limit = CAP - ItemSize[JStar];
  for (i=1; i<JStar; i++)
  {
    if (ItemSize[i] <= Limit)
    {
      JPrime = i;
      break;
    }
  }

  if (JPrime == 0) JPrime = JStar;

  CJ2 = JStar - JPrime;

  SJ2 = 0;
  for (i=JPrime; i<JStar; i++) SJ2 += ItemSize[i];

  JDPrime = JStar;
  SJ3 = ItemSize[JDPrime];

  if (JDPrime < n)
  {
    while (ItemSize[JDPrime+1] == ItemSize[JDPrime])
    {
      JDPrime++;
      SJ3 += ItemSize[JDPrime];
      if (JDPrime == n) break;
    }
  }

  L2 = CJ12;

  do
  {
    AddSizeSum = SJ3 + SJ2 - (CJ2 * CAP);
    AddBins = 0;
    AddCAP = 0;

    while (AddCAP < AddSizeSum)
    {
      AddBins++;
      AddCAP += CAP;
    }

    if ((CJ12 + AddBins) > L2)
    L2 = (CJ12 + AddBins);

    JDPrime++;
    if (JDPrime > n) break;

    SJ3 += ItemSize[JDPrime];

    if (JDPrime < n)
    {
      while (ItemSize[JDPrime+1] == ItemSize[JDPrime])
      {
        JDPrime++;
        SJ3 += ItemSize[JDPrime];
        if (JDPrime == n) break;
      }
    }

    if (JPrime > 1)
    {
      Limit = CAP - ItemSize[JDPrime];
      while (ItemSize[JPrime-1] <= Limit)
      {
        JPrime--;
        CJ2++;
        SJ2 += ItemSize[JPrime];
        if (JPrime == 1) break;
      }
    }
  } while (JDPrime <= n);

  *LB = L2;

  /* printf("BP_MTL2: L2 = %d\n",L2); */
  /* printf("BP_MTL2: END\n"); */
}

void BP_ModifiedMTL3(double CAP,
                     double *ItemSize,
                     int n,
                     int UB,
                     int *LB)
{
  int i,j,k,nBar;
  int L3,Z,Zr,L2;
  int *Bin;
  double *W;

  /*
  printf("BP_ModifiedMTL3: Input:\n");
  printf("UB=%d, CAP=%d, n=%d\n",UB,CAP,n);
  printf("ItemSize =");
  for (i=1; i<=n; i++) printf(" %d",ItemSize[i]);
  printf("\n");
  */

  W = MemGetDV(n+1);
  Bin = MemGetIV(n+1);

  for (i=1; i<=n; i++) W[i] = ItemSize[i];

  L3 = 0;
  Z = 0;
  nBar = n;

  while (nBar >= 1)
  {
    BP_DominancePacking(CAP,W,nBar,Bin,&Zr);

    if (Zr > 0)
    {
      Z += Zr;

      k = 0;
      for (j=1; j<=nBar; j++)
      {
        if (Bin[j] == 0)
        {
          k++;
          W[k] = W[j];
        }
      }

      nBar = k;
    }

    if (nBar == 0)
    L2 = 0;
    else
    BP_MTL2(CAP,W,nBar,&L2);

    if (L3 < (Z + L2))
    L3 = (Z + L2);

    if (L3 >= UB) break;

    /* nBar--; */
    nBar = 0; /* Skip the relaxations */
  }

  *LB = L3;

  MemFree(W);
  MemFree(Bin);

  /* printf("END BP_ModifiedMTL3 = %d\n",L3); */
}

//NOTE: This assumes that ItemSize is sorted decreasing.
void BP_ExactBinPacking(double CAP,
                        double *ItemSize,
                        int n,
                        int *LB,
                        int *UB,
                        int *Bin)
{
   char ComputeBound;
   int i,j,Item,ThisBin,NewBin,LoopNr;
   double TmpMaxSize,AccSlack,MaxTotalSlack,TotalItemSizeSum,MinItemSize;
   int FirstCandidateBin;
   int TmpIndex;
   int LowerBoundBins,UpperBoundBins,TmpLowerBound,TmpUpperBound;
   int CompressedItems,CompressedBins;
   int *BestBin, *B, *TmpBin;
   int **MinItem;
   double *W,*BinSlack;


   //printf("BEGIN BP_ExactBinPacking(NoOfItems=%d, CAP=%d)\n",n,CAP);
   //printf("ItemSize =");
   //for (i=1; i<=n; i++) printf(" %d",ItemSize[i]);
   //printf("\n");
  

   BestBin = MemGetIV(n+1);
   B = MemGetIV(n+1);
   BinSlack = MemGetDV(n+1);
   W = MemGetDV(n+1);
   TmpBin = MemGetIV(n+1);

   BP_FirstFit(CAP,ItemSize,n,BestBin,&UpperBoundBins);
   for (i=1; i<=n; i++) Bin[i] = BestBin[i];

   //printf("UpperBoundBins = %d\n",UpperBoundBins);

   MinItem = MemGetIM(n+1,UpperBoundBins+1);

   //printf("Calling MTL2\n"); 
   BP_MTL2(CAP,ItemSize,n,&LowerBoundBins);

   //printf("MTL2 = %d\n",LowerBoundBins);

   /*
     if (LowerBoundBins < UpperBoundBins)
     BP_ModifiedMTL3(CAP,ItemSize,n,UpperBoundBins,&LowerBoundBins);
   */

   LoopNr = 0;

   if (LowerBoundBins == UpperBoundBins)
      {
         goto EndOfExactBinPacking;
      }

   /* If we get to here, UpperBoundBins is at least 3. */

   MinItemSize = ItemSize[n];
   /* printf("MinItemSize = %d\n",MinItemSize); */
   //if (MinItemSize <= 0) exit(0);

   TotalItemSizeSum = 0;
   for (i=1; i<=n; i++) TotalItemSizeSum += ItemSize[i];

   /* printf("TotalItemSizeSum = %d\n",TotalItemSizeSum); */

   MaxTotalSlack = ((UpperBoundBins-1) * CAP) - TotalItemSizeSum;
   /* Max. slack, if a better solution than the currently best solution
      is to be found. */

   /* printf("MaxTotalSlack = %d\n",MaxTotalSlack); */

   /* Find the minimum level at which a different allocation must be
      made if a better solution than the currently best exists. */

   for (i=1; i<=n; i++) BinSlack[i] = CAP;

   AccSlack = 0;
   Item = 0;

   for (i=1; i<=n-1; i++)
      {
         ThisBin = BestBin[i];
         B[i] = ThisBin;
         BinSlack[ThisBin] -= ItemSize[i];

         if (BinSlack[ThisBin] < MinItemSize)
            {
               AccSlack += BinSlack[ThisBin];
               if (AccSlack > MaxTotalSlack)
                  {
                     Item = i;
                     break;
                  }
            }

         if (ThisBin <= UpperBoundBins-2)
            {
               Item = i;
            }
      }

   /* printf("Initial Item = %d, AccSlack=%d\n",Item,AccSlack); */

   /* Redo the stuff until this Item */
   for (i=1; i<=n; i++) BinSlack[i] = CAP;
   for (i=1; i<=n; i++) B[i] = 0;

   AccSlack = 0;

   for (i=1; i<=Item; i++)
      {
         ThisBin = BestBin[i];
         B[i] = ThisBin;
         /* printf("B(%d):= %d\n",i,B[i]); */
         BinSlack[ThisBin] -= ItemSize[i];
         /* printf("BinSlack(%d):= %d\n",ThisBin,BinSlack[ThisBin]); */
         if (BinSlack[ThisBin] < MinItemSize)
            AccSlack += BinSlack[ThisBin];
      }

   /* printf("Item = %d, AccSlack=%d\n",Item,AccSlack); */

   for (i=0; i<=n; i++)
      for (j=1; j<=UpperBoundBins; j++)
         MinItem[i][j] = i+1;

   BinSlack[0] = -1; /* => Different from CAP */

   LoopNr = 0;

   do
      {
         LoopNr++;
         /* Reallocate item */
         //printf("Top of do loop...: LoopNr=%d Item=%d\n",LoopNr,Item); 
         /* Remove the item from its current bin */
         ThisBin = B[Item];

         /* printf("Item=%d, ThisBin=%d\n",Item,ThisBin); */

         if (ThisBin > 0)
            if (BinSlack[ThisBin] < MinItemSize)
               {
                  /*
                    printf("Item closed bin: BinSlack(%d)=%d, MinItemSize=%d\n",
                    ThisBin,BinSlack[ThisBin],MinItemSize);
                  */

                  /* This item closed ThisBin when it was allocated to this bin */
                  /* Some dominated allocations may be derived */

                  AccSlack -= BinSlack[ThisBin];

                  /* printf("AccSlack:= %d\n",AccSlack); */

                  /* Find the first item j for which w(j)+w(n) fits into the bin
                     instead of the current item. Item j is the first nondominated
                     item that can be filled into this bin, given the allocation
                     of items 1,...,Item-1. */

                  TmpMaxSize = BinSlack[ThisBin] + ItemSize[Item] - ItemSize[n];
                  TmpIndex = n+1;

                  /* printf("TmpMaxSize=%d\n",TmpMaxSize); */

                  for (j=Item+1; j<=n; j++)
                     if (ItemSize[j] <= TmpMaxSize)
                        {
                           TmpIndex = j;
                           break;
                        }

                  MinItem[Item-1][ThisBin] = TmpIndex;
                  /* printf("MinItem(%d,%d):= %d\n",Item-1,ThisBin,MinItem[Item-1][ThisBin]); */
               }

         if (ThisBin > 0)
            BinSlack[ThisBin] += ItemSize[Item];

         /* printf("BinSlack(%d):= %d\n",ThisBin,BinSlack[ThisBin]); */

         /* Find the next bin > B[Item] for this item */
         if (BinSlack[ThisBin] == CAP)
            {
               /* The item initialized a bin, so it should not be tried in
                  any other bin: move up the tree to the previous item */

               /* printf("Item initialized new bin => NewBin:= 0\n"); */
               NewBin = 0;
               /* Backtrack */
            }
         else
            {
               /* Find the first bin > B[Item] for which the slack before
                  allocation of this item into the bin differs from *all*
                  smaller-indexed bins. This is also a dominance criterion.
                  It may be a new criterion (!).
                  The criterion follows from the observation that there is
                  no point in assigning this item into a bin which effectively
                  is a copy of a smaller-indexed bin.

                  The new bin for this item must also satisfy the dominance
                  criterion wrt. the MinItem.

                  In addition, if w(j) = w(j+1), then item j+1 is required to
                  satisfy B(j+1) >= B(j). This is also a new criterion.
               */

               FirstCandidateBin = ThisBin+1;

               /* printf("FirstCandidateBin:= %d\n",FirstCandidateBin); */

               if ((ItemSize[Item] == ItemSize[Item-1]) &&
                   (FirstCandidateBin < B[Item-1]))
                  {
                     FirstCandidateBin = B[Item-1];
                     /*
                       printf("Identical items %d and %d: Adj. FirstCandidateBin:= %d\n",
                       Item,Item-1,FirstCandidateBin);
                     */
                  }

               NewBin = 0;

               for (j=FirstCandidateBin; j<=UpperBoundBins-1; j++)
                  {
                     NewBin = j;
                     /* printf("Possible bin %d:\n",j); */

                     /* printf("MinItem(%d,%d) = %d\n",Item-1,j,MinItem[Item-1][j]); */

                     if (BinSlack[j] < ItemSize[Item])
                        {
                           /*
                             printf("ItemSize=%d > BinSlack(%d)=%d\n",
                             ItemSize[Item],j,BinSlack[j]);
                           */
                           NewBin = 0;
                           continue;
                        }

                     if (Item < MinItem[Item-1][j])
                        {
                           /* Try next bin */
                           /*
                             printf("Item < MinItem(%d,%d) = %d\n",Item-1,j,MinItem[Item-1][j]);
                           */
                           NewBin = 0;
                           continue;
                        }

                     if (NewBin > 0)
                        {
                           /* printf("Check acc. slack:\n"); */
                           /* Check the accumulated BinSlack */
                           if ((BinSlack[NewBin] - ItemSize[Item]) < MinItemSize)
                              {
                                 /*
                                   printf("Item generates slack of %d in bin %d\n",
                                   BinSlack[NewBin] - ItemSize[Item],NewBin);
                                 */
                                 /* The Item will close the NewBin */
                                 if ((AccSlack +
                                      (BinSlack[NewBin] - ItemSize[Item])) > MaxTotalSlack)
                                    {
                                       /* printf(" + AccSlack = %d > %d\n",AccSlack,MaxTotalSlack); */
                                       NewBin = 0;
                                    }
                              }
                        }

                     if (NewBin > 0)
                        {
                           for (i=1; i<j; i++)
                              {
                                 if (BinSlack[i] == BinSlack[j])
                                    {
                                       /* Try next bin */

                                       /*
                                         printf("Equal bin slacks (tried bin %d):\n",j);
                                         printf("BinSlack(%d) = BinSlack(%d) = %d\n",
                                         i,j,BinSlack[i]);
                                       */

                                       NewBin = 0;
                                       break;
                                    }
                              }
                        }

                     if (NewBin > 0) break;
                  }
            }

         //printf("Item=%d, NewBin:= %d\n",Item,NewBin);
         
         if (NewBin > 0)
            {
               /* Assign Item to NewBin */

               B[Item] = NewBin;
               BinSlack[NewBin] -= ItemSize[Item];

               /*
                 printf("B(%d):= %d, BinSlack(%d):=%d\n",
                 Item,B[Item],NewBin,BinSlack[NewBin]);
               */

               if (BinSlack[NewBin] < MinItemSize)
                  {
                     AccSlack += BinSlack[NewBin];
                     /* printf("AccSlack:= %d\n",AccSlack); */
                  }

               /* Compute new lower bound given this allocation */
               ComputeBound = 0;

               CompressedItems = 0;
               for (i=1; i<=UpperBoundBins; i++)
                  {
                     if (BinSlack[i] < CAP)
                        {
                           W[++CompressedItems] = (CAP - BinSlack[i]);
                           /*
                             printf("From Bins: W(%d):= %d\n",
                             CompressedItems,W[CompressedItems]);
                           */
                           if (W[CompressedItems] != ItemSize[CompressedItems])
                              ComputeBound = 1;
                        }
                  }

               CompressedBins = CompressedItems;

               for (i=Item+1; i<=n; i++)
                  {
                     W[++CompressedItems] = ItemSize[i];
                     /*
                       printf("From Item %d: W(%d):= %d\n",
                       i,CompressedItems,W[CompressedItems]);
                     */
                     /* TmpBin[CompressedItems] = 0; */
                  }

               /* Compute new feasible solution given this allocation */
               for (i=1; i<=CompressedBins; i++) TmpBin[i] = i;
               for (i=CompressedBins+1; i<=CompressedItems; i++) TmpBin[i] = 0;

               BP_ModifiedFirstFit(CAP,W,CompressedItems,TmpBin,&TmpUpperBound);
               /* printf("FFD => %d bins\n",TmpUpperBound); */

               if (TmpUpperBound < UpperBoundBins)
                  {
                     UpperBoundBins = TmpUpperBound;
                     /* printf("UpperBoundBins:= %d\n",UpperBoundBins); */

                     /* for (i=1; i<=Item; i++) printf("B(%d) = %d\n",i,B[i]); */

                     /* for (i=1; i<=n; i++) Bin[i] = TmpBin[i]; */
                     for (i=Item+1; i<=n; i++) Bin[i] = TmpBin[i-Item+CompressedBins];
                     for (i=1; i<=Item; i++) Bin[i] = B[i];

                     if (LowerBoundBins == UpperBoundBins)
                        goto EndOfExactBinPacking;
                  }

               ComputeBound = 1;

               /* Sort items from bins */
               /* SortIVDec(W,CompressedItems); */
               if (CompressedBins > 1)
                  SortDVDec(W,CompressedBins);
               /* Sufficient to sort the compressed bins */

               if (ComputeBound == 0)
                  {
                     /* printf("ComputeBound = 0\n"); */
                     TmpLowerBound = LowerBoundBins;
                     /* printf("TmpLowerBound = %d\n",TmpLowerBound); */
                  }
               else
                  {
                     BP_MTL2(CAP,W,CompressedItems,&TmpLowerBound);

                     if (TmpLowerBound < UpperBoundBins)
                        BP_ModifiedMTL3(CAP,W,CompressedItems,UpperBoundBins,&TmpLowerBound);

                     /* printf("TmpLowerBound = %d\n",TmpLowerBound); */
                  }


               if (TmpLowerBound >= UpperBoundBins)
                  {
                     /* Try next bin for this item in the next loop */
                     /* No action is OK: keep the current item for reallocation to the
                        next bin */
                  }
               else
                  {
                     /* Forward step = branch */

                     for (i=1; i<=UpperBoundBins; i++)
                        MinItem[Item][i] = MinItem[Item-1][i];

                     B[Item+1] = 0;
                     Item++;
                  }
            }
         else
            {
               /* NewBin = 0 */
               /* No more bins need to be tried for this item */
               /* Backtrack */
               //printf("Item:%d set to 0\n",Item);
               B[Item] = 0;
               Item--;
            }

         /* if ((LoopNr % 10) == 0) printf("LoopNr = %d\n",LoopNr); */

      } while ((LoopNr < 1000) && (Item > 1));

   if (Item == 1)
      LowerBoundBins = UpperBoundBins;

 EndOfExactBinPacking:

   *LB = LowerBoundBins;
   *UB = UpperBoundBins;

   /*
     printf("LoopNr=%d, LowerBoundBins=%d, UpperBoundBins=%d\n",
     LoopNr,LowerBoundBins,UpperBoundBins);
   */

   MemFree(BestBin);
   MemFree(B);
   MemFree(BinSlack);
   MemFree(W);
   MemFree(TmpBin);

   MemFreeIM(MinItem,n+1);

   return;
}

