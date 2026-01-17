/*
 * Copyright (c) 2025 Yu Yang & Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */


#include <algorithm>
#include <numeric>
#include "route_opt_macro.hpp"
#include "deluxing_macro.hpp"
#include "deluxing_controller.hpp"
#include "deluxing_kmeans.hpp"

namespace RouteOpt::DeLuxing {
    void DeLuxingController::deLuxing(Solver &solver, double UB, int NClust, int beta1, int beta2,
                                      std::vector<int> &Idxdel,
                                      double Timelimit, double Tolerance, bool Verbose) {
        auto beg = std::chrono::high_resolution_clock::now();

        // record the original solver environment
        int method, crossover, threads;
        SAFE_SOLVER(solver.getEnvMethod(&method))
        SAFE_SOLVER(solver.getCrossOver(&crossover))
        SAFE_SOLVER(solver.getThreads(&threads))

        // set the solver environment for optimization
        SAFE_SOLVER(solver.setEnvMethod(SOLVER_BARRIER))
        SAFE_SOLVER(solver.setEnvCrossOver(SOLVER_CROSSOVER_DOWN))
        SAFE_SOLVER(solver.setEnvThreads(NUM_THREADS, false))
        SAFE_SOLVER(solver.reoptimize())

        int NRow, NCol;
        SAFE_SOLVER(solver.getNumRow(&NRow))
        SAFE_SOLVER(solver.getNumCol(&NCol))
        int NN = NCol;

        size_t numnzP;
        SAFE_SOLVER(solver.XgetVars(&numnzP, nullptr, nullptr, nullptr, 0, NCol))
        std::vector<double> dj(NCol);


        int nv = std::max(NCol, NRow);
        numnzP = (numnzP < nv ? nv : numnzP);

        // For double type vectors
        std::vector<double> cobj(nv),
                obj_real(NCol),
                cmatval(numnzP),
                matval(numnzP),
                rhs_real(NRow),
                dual(NRow);

        // For size_t type vectors
        std::vector<size_t> cmatbeg(nv + 1),
                matbeg(NCol + 1);

        // For int type vectors
        std::vector<int> cmatind(numnzP),
                matind(numnzP),
                del_idx(NCol);

        // For bool type vectors
        std::vector<bool> ifdel(NCol),
                delall(NCol);

        SAFE_SOLVER(solver.getRhs(0, NRow, rhs_real.data()))
        SAFE_SOLVER(solver.XgetVars(&numnzP, matbeg.data(), matind.data(), matval.data(), 0, NCol))
        matbeg[NCol] = numnzP;


        double LB;
        SAFE_SOLVER(solver.getObjVal(&LB))
        double gap = UB - LB + Tolerance, gap_half = (UB - LB) / 2.0 + Tolerance;
        int num_half = 0, ncol_del0 = 0;
        SAFE_SOLVER(solver.XgetVars(&numnzP, cmatbeg.data(), cmatind.data(), cmatval.data(), 0, NCol))
        matbeg[NCol] = cmatbeg[NCol] = numnzP;
        SAFE_SOLVER(solver.getRC(0, NCol, dj.data()))
        SAFE_SOLVER(solver.getObj(0, NCol, cobj.data()))
        std::vector<std::pair<double, int> > RCVec(NCol);
        for (int k = 0; k < NCol; ++k)
            RCVec[k] = std::make_pair(dj[k], k);
        std::stable_sort(RCVec.begin(), RCVec.end());


        std::vector<int> GlobalPool(NCol);
        std::iota(GlobalPool.begin(), GlobalPool.end(), 0);

        std::vector<int> PoolCopy = GlobalPool;
        for (int k = 0; k < NCol; ++k) {
            GlobalPool[k] = PoolCopy[RCVec[k].second];
            if (dj[k] > gap) ++ncol_del0;
        }

        matbeg[0] = numnzP = 0;
        for (int k = 0; k < NCol - ncol_del0; ++k) {
            if (RCVec[k].first > gap_half && !num_half)
                num_half = k;
            for (auto j = cmatbeg[RCVec[k].second]; j < cmatbeg[RCVec[k].second + 1]; ++j) {
                matind[numnzP] = cmatind[j];
                matval[numnzP++] = cmatval[j];
            }
            matbeg[k + 1] = numnzP;
            obj_real[k] = cobj[RCVec[k].second];
        }
        std::vector<int> idx(NCol);
        std::vector<uint32_t> clusters(NCol, 0);
        if (ncol_del0 >= NCol) goto Terminate;
        std::cout << "|R| = " << NN << ", initial del = " << ncol_del0 << ", |R_1| = " << num_half << std::endl;
        std::cout << "------------------------------------------------------------------------------------" <<
                std::endl;


        std::iota(idx.data(), idx.data() + NCol, 0);
        SAFE_SOLVER(solver.delVars(NCol, idx.data()))
        SAFE_SOLVER(solver.updateModel())

        GlobalPool.erase(GlobalPool.begin() + NCol - ncol_del0, GlobalPool.end());
        NCol -= ncol_del0;
        GlobalPool.resize(NCol);
        PoolCopy.clear();
        RCVec.clear();
        SAFE_SOLVER(solver.XaddVars(num_half, matbeg[num_half], matbeg.data(), matind.data(), matval.data(),
            obj_real.data(), nullptr, nullptr, nullptr, nullptr))

        cmatbeg[0] = 0;
        for (int k = 0; k < NRow; ++k) {
            cmatind[k] = k;
            cmatval[k] = -1.;
            cobj[k] = UB;
            cmatbeg[k + 1] = k + 1;
        }
        SAFE_SOLVER(solver.XaddVars(NRow, NRow, cmatbeg.data(), cmatind.data(), cmatval.data(), cobj.data(),
            nullptr, nullptr, nullptr, nullptr))

        if (NCol <= beta1) {
            std::vector<std::vector<double> > data;
            for (int j = 0; j < NCol; ++j) {
                std::vector<double> col(rhs_real.data(), rhs_real.data() + NRow);
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    col[matind[k]] -= matval[k];
                data.emplace_back(col);
            }
            KMeans::clustering_parameters<double> parameters(NClust);
            parameters.set_random_seed(0);
            parameters.set_max_iteration(100);
            std::vector<std::vector<double> > means;
            kmeans_lloyd_parallel(data, means, clusters, parameters);
        } else {
            std::vector<std::pair<double, int> > col_norm(NCol);
            for (int j = 0; j < NCol; ++j) {
                for (int k = 0; k < NRow; ++k)
                    cmatval[k] = rhs_real[k];
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    if (matind[k] < NRow)
                        cmatval[matind[k]] -= matval[k];
                double sum = 0.;
                for (int k = 0; k < NRow; ++k)
                    sum += cmatval[k] * cmatval[k];
                col_norm[j] = std::make_pair(sum, j);
            }
            std::stable_sort(col_norm.begin(), col_norm.end());
            double rate = 1.0 / NClust;
            for (int r = 0; r < NClust; ++r) {
                int s1 = static_cast<int>(r * NCol * rate), s2 = (r == NClust - 1
                                                                      ? NCol
                                                                      : static_cast<int>((r + 1) * NCol * rate));
                for (int j = s1; j < s2; ++j)
                    clusters[col_norm[j].second] = r;
            }
        }

        for (int r = 0; r < NClust; ++r) {
            int ncol_del = 0;
            for (int k = 0; k < NRow; ++k) cmatval[k] = 0.;

            for (int j = 0; j < NCol; ++j)
                if (clusters[j] == r) {
                    for (int k = 0; k < NRow; ++k)
                        cmatval[k] += rhs_real[k];
                    for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                        cmatval[matind[k]] -= matval[k];
                }


            SAFE_SOLVER(solver.setRhs(0, NRow, cmatval.data()))
            SAFE_SOLVER(solver.reoptimize())
            SAFE_SOLVER(solver.getDual(0, NRow, dual.data()))

            double threshold = UB + Tolerance;
            for (int k = 0; k < NRow; ++k)
                threshold -= dual[k] * rhs_real[k];

            double minrc = 0.;
            for (int j = num_half; j < NCol; ++j) {
                double rcj = obj_real[j];
                for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                    rcj -= dual[matind[k]] * matval[k];
                if (rcj >= threshold) {
                    ifdel[j] = true;
                    ++ncol_del;
                } else {
                    ifdel[j] = false;
                    if (rcj < minrc) minrc = rcj;
                }
            }
            threshold -= minrc;
            SAFE_SOLVER(solver.getRC(0, num_half, dj.data()))
            int dc = 0, real_half = num_half;
            std::vector<int> cp_(real_half);
            for (int j = 0; j < num_half; ++j) {
                cp_[j] = j;
                if (dj[j] >= threshold) {
                    ifdel[j] = true;
                    ++ncol_del;
                    del_idx[dc++] = j;
                } else ifdel[j] = false;
            }


            if (dc) {
                for (int j = 0; j < num_half; ++j)
                    if (ifdel[j])
                        cp_.erase(std::remove(cp_.begin(), cp_.end(), j), cp_.end());
                real_half -= dc;
                SAFE_SOLVER(solver.delVars(dc, del_idx.data()))
                SAFE_SOLVER(solver.updateModel())
            }

            for (int j = 0; j < NCol; ++j) {
                delall[j] = ifdel[j];
                ifdel[j] = ifdel[j] && clusters[j] != r;
            }
            int extra = ncol_del, ite = 1;
            while (extra > beta2) {
                extra = 0;
                for (int k = 0; k < NRow; ++k) cmatval[k] = 0.;
                for (int j = 0; j < NCol; ++j) {
                    if (ifdel[j]) {
                        for (int k = 0; k < NRow; ++k)
                            cmatval[k] += rhs_real[k];
                        for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                            cmatval[matind[k]] -= matval[k];
                    }
                }
                SAFE_SOLVER(solver.setRhs(0, NRow, cmatval.data()))
                SAFE_SOLVER(solver.reoptimize())
                SAFE_SOLVER(solver.getDual(0, NRow, dual.data()))

                threshold = UB + Tolerance;
                for (int k = 0; k < NRow; ++k)
                    threshold -= dual[k] * rhs_real[k];

                minrc = 0.;
                for (int j = num_half; j < NCol; ++j) {
                    ifdel[j] = false;
                    if (delall[j]) continue;
                    double rcj = obj_real[j];
                    for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                        rcj -= dual[matind[k]] * matval[k];
                    if (rcj >= threshold) {
                        delall[j] = ifdel[j] = true;
                        ++extra;
                    } else if (rcj < minrc) minrc = rcj;
                }
                threshold -= minrc;
                for (int j = 0; j < num_half; ++j) {
                    ifdel[j] = false;
                    if (delall[j]) continue;
                    double rcj = obj_real[j];
                    for (auto k = matbeg[j]; k < matbeg[j + 1]; ++k)
                        rcj -= dual[matind[k]] * matval[k];
                    if (rcj >= threshold) {
                        delall[j] = ifdel[j] = true;
                        ++extra;
                    }
                }

                dc = 0;
                for (int j = 0; j < num_half; ++j)
                    if (ifdel[j])
                        del_idx[dc++] = static_cast<int>(std::find(cp_.begin(), cp_.end(), j) - cp_.begin());

                if (dc) {
                    for (int j = 0; j < num_half; ++j)
                        if (ifdel[j])
                            cp_.erase(std::remove(cp_.begin(), cp_.end(), j), cp_.end());
                    real_half -= dc;
                    SAFE_SOLVER(solver.delVars(dc, del_idx.data()))
                    SAFE_SOLVER(solver.updateModel())
                }

                ncol_del += extra;
                ++ite;
                if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - beg).count() > Timelimit)
                    break;
            }


            int kc = 0;
            numnzP = 0;
            for (int i = 0; i < NCol; ++i) {
                if (!delall[i]) {
                    for (auto k = matbeg[i]; k < matbeg[i + 1]; ++k) {
                        matind[numnzP] = matind[k];
                        matval[numnzP++] = matval[k];
                    }
                    obj_real[kc] = obj_real[i];
                    GlobalPool[kc] = GlobalPool[i];
                    clusters[kc] = clusters[i];
                    matbeg[++kc] = numnzP;
                }
            }

            num_half = real_half;
            GlobalPool.resize(kc);
            NCol -= ncol_del;
            ncol_del += ncol_del0;

            if (Verbose)
                std::cout << "Round " << r + 1 << ", deleted " << ncol_del << " cols, ite = " << ite << ",  #cols = " <<
                        NCol << std::endl;
            ncol_del0 = 0;
        }

        SAFE_SOLVER(solver.getNumCol(&NCol))
        std::iota(idx.data(), idx.data() + NCol, 0);
        SAFE_SOLVER(solver.delVars(NCol, idx.data()))
        SAFE_SOLVER(solver.updateModel())
        NCol = static_cast<int>(GlobalPool.size());
        SAFE_SOLVER(solver.XaddVars(NCol, numnzP, matbeg.data(), matind.data(), matval.data(),
            obj_real.data(), nullptr, nullptr, nullptr, nullptr))
        SAFE_SOLVER(solver.setRhs(0, NRow, rhs_real.data()))
        SAFE_SOLVER(solver.updateModel())

    Terminate:
        // recover the original solver environment
        SAFE_SOLVER(solver.setEnvMethod(method))
        SAFE_SOLVER(solver.setEnvCrossOver(crossover))
        SAFE_SOLVER(solver.setEnvThreads(threads, false))

        std::stable_sort(GlobalPool.begin(), GlobalPool.end());

        if (NCol) {
            for (int j = 0; j < GlobalPool[0]; ++j)
                Idxdel.push_back(j);
            for (int k = 0; k < NCol - 1; ++k)
                for (int j = GlobalPool[k] + 1; j < GlobalPool[k + 1]; ++j)
                    Idxdel.push_back(j);
            for (int j = GlobalPool[NCol - 1] + 1; j < NN; ++j)
                Idxdel.push_back(j);
        }
        auto eps = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - beg).count();
        std::cout << SMALL_PHASE_SEPARATION;
        printTimeMessage("deluxing", eps);
        std::cout << "|R| = " << NN << ", NRemain = " << NCol <<
                ", percentage = "
                << (
                    NN - NCol) * 100.0 / NN << "%\n\n" << std::endl;
    }
}
