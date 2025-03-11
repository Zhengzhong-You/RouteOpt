/* 
 * Copyright (c) 2025 Zhengzhong (Ricky) You.
 * All rights reserved.
 * Software: RouteOpt
 * License: GPL-3.0
 */

#include "cvrp_pricing_controller.hpp"

namespace RouteOpt::Application::CVRP {
    void CVRP_Pricing::freeMemory() const {
        delete[]all_label;
        delete[]col_pool4_pricing;
        delete[]copy_col_pool4_pricing;

        if (rc2_till_this_bin_in_forward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]rc2_till_this_bin_in_forward_sense[i];
            }
            delete[]rc2_till_this_bin_in_forward_sense;
        }

        if (rc2_bin_in_forward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]rc2_bin_in_forward_sense[i];
            }
            delete[]rc2_bin_in_forward_sense;
        }

        if (label_array_in_forward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]label_array_in_forward_sense[i];
            }
            delete[]label_array_in_forward_sense;
        }


        if (if_exist_extra_labels_in_forward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]if_exist_extra_labels_in_forward_sense[i];
            }
            delete[]if_exist_extra_labels_in_forward_sense;
        }

        //backwards
        if (rc2_till_this_bin_in_backward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]rc2_till_this_bin_in_backward_sense[i];
            }
            delete[]rc2_till_this_bin_in_backward_sense;
        }
        if (rc2_bin_in_backward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]rc2_bin_in_backward_sense[i];
            }
            delete[]rc2_bin_in_backward_sense;
        }
        if (label_array_in_backward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]label_array_in_backward_sense[i];
            }
            delete[]label_array_in_backward_sense;
        }
        if (if_exist_extra_labels_in_backward_sense != nullptr) {
            for (int i = 0; i < dim; ++i) {
                delete[]if_exist_extra_labels_in_backward_sense[i];
            }
            delete[]if_exist_extra_labels_in_backward_sense;
        }
    }
}
