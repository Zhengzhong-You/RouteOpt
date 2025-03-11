/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_FCITS
#define _H_FCITS

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
                      CnstrMgrPointer CutsCMP);

#endif

