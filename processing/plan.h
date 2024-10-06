#ifndef PROCESSING_PLAN_H
#define PROCESSING_PLAN_H

#include <string>

#include "utils/config.h"
#include "graph/graph.h"
#include "utils/nucleus/nd_interface.h"

struct InitialOrder
{
    uint8_t u[MAX_VCOUNT];
};

struct Plan
{
    uint8_t start_vertex_id_;
    uint8_t mask_size_; // the total number of groups
    // each mask[i] is a group of vertices. mask[0] is the start vertex, mask[1] is the first group of vertices, etc.
    // If vertices v in the group of mask[i], the v bit in mask[i] is set to 1.
    uint16_t masks_[MAX_VCOUNT];
    // mask_size_prefix_sum_[i] is the total number of vertices in the first i-1 groups. (not include the ith group)
    uint8_t mask_size_prefix_sum_[MAX_VCOUNT];
    uint8_t res_pos[MAX_VCOUNT];  // order the vertex according to the group and query id. if mask[0]="10", then res_pos[2]=0
    uint16_t is_tree_group_; // if ith group is a tree group, then the ith bit in is_tree_group_ is set to 1.

    Plan() {}

    Plan(
        const Graph& query,
        uint32_t *relation_sizes,
        std::string method
    );

    Plan(
        const Graph& query,
        uint32_t *relation_sizes,
        uint32_t gsi
    );

    void AddGroup(uint8_t i, InitialOrder& initial_order);

    void Print(const Graph& query);
};

extern __constant__ Plan c_plan;


#endif //PROCESSING_PLAN_H
