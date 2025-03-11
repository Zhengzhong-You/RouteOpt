/* SAS modified this file. */
/* (C) Copyright 2003 Jens Lysgaard. All rights reserved. */
/* OSI Certified Open Source Software */
/* This software is licensed under the Common Public License Version 1.0 */

#ifndef _H_NEWHTOUR
#define _H_NEWHTOUR

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
                            double *MaxHTourViolation);

#endif

