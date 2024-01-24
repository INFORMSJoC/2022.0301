// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef MIP_GRAPH_HPP
#define MIP_GRAPH_HPP

#include <gurobi_c++.h>                     // for GRBVar
#include <boost/graph/adjacency_list.hpp>   // for source, vecS (ptr only)
#include <boost/graph/graph_selectors.hpp>  // for bidirectionalS
#include <boost/graph/graph_traits.hpp>     // for graph_traits, graph_trait...
#include <boost/pending/property.hpp>       // for lookup_one_property_inter...
#include <cstddef>                          // for size_t
#include <ostream>                          // for operator<<, ostream, basi...
#include "Job.h"                            // for Job
#include "ModernDD/NodeBddTable.hpp"        // for NodeTableEntity
#include "ModernDD/NodeId.hpp"              // for NodeId
#include "NodeBdd.hpp"                      // for NodeBdd, NodeBdd::dbl_array

struct VertexData {
    size_t index{};
    NodeId node_id{};
};

struct EdgeData {
    size_t id{};
    bool   high{};
    GRBVar x{};
};
using MipGraph = boost::adjacency_list<boost::vecS,
                                       boost::vecS,
                                       boost::bidirectionalS,
                                       VertexData,
                                       EdgeData>;
using Edge = boost::graph_traits<MipGraph>::edge_descriptor;
using Vertex = boost::graph_traits<MipGraph>::vertex_descriptor;

class ColorWriterEdgeX {
   private:
    MipGraph&                 g;
    NodeTableEntity<NodeBdd>* table;
    static constexpr double   EPS_GRAPH = 1e-6;

   public:
    explicit ColorWriterEdgeX(MipGraph& _g, NodeTableEntity<NodeBdd>* _table)
        : g{_g},
          table(_table) {}

    void operator()(std::ostream& output, const Edge& _edge) {
        auto  node_id = g[source(_edge, g)].node_id;
        auto& node = table->node(node_id);

        if (g[_edge].high) {
            auto& x = node.get_lp_x()[1];
            if (x > EPS_GRAPH) {
                output << "[label = " << x << ",color = red]";
            } else {
                output << "[label = " << x << "]";
            }
        } else {
            auto& x = node.get_lp_x()[0];
            if (x > EPS_GRAPH) {
                output << "[label = " << x << ",color = red, style = dashed]";
            } else {
                output << "[label = " << x << ",style=dashed]";
            }
        }
    }
};

class ColorWriterEdgeIndex {
   private:
    const MipGraph& g;

   public:
    explicit ColorWriterEdgeIndex(const MipGraph& _g) : g{_g} {}

    void operator()(std::ostream& output, const Edge& _edge) {
        auto index = g[_edge].id;
        auto high = g[_edge].high;

        if (high) {
            output << "[label = " << index << "]";
        } else {
            output << "[label = " << index << ", style = dashed]";
        }
    }
};

class ColorWriterVertex {
   private:
    const MipGraph&                 g;
    const NodeTableEntity<NodeBdd>& table;

   public:
    ColorWriterVertex(const MipGraph&                 _g,
                      const NodeTableEntity<NodeBdd>& _table)
        : g{_g},
          table{_table} {}

    void operator()(std::ostream& output, const Vertex& _vertex) {
        const auto& tmp_node_id = g[_vertex].node_id;
        if (tmp_node_id > 1) {
            output << " [label=\" " << table.node(tmp_node_id).get_job()->job
                   << " " << table.node(tmp_node_id).get_weight() << "\"]";
        }
    }
};

#endif  // MIP_GRAPH_HPP
