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

#ifndef LABEL_HPP
#define LABEL_HPP

#include <limits>               // for numeric_limits
#include "ModernDD/NodeId.hpp"  // for NodeId
struct Job;

template <typename N>
class Label {
   private:
    int       weight{-1};
    bool      high{false};
    Label<N>* prev_label{nullptr};
    NodeId    node_id{};

    double f{std::numeric_limits<double>::max()};
    Job*   label_job{nullptr};
    Job*   prev_job{nullptr};

   public:
    /**
     * Constructor
     */
    Label() = default;
    Label(const Label<N>& src) = default;
    Label(Label<N>&& src) noexcept = default;
    auto operator=(const Label<N>& src) -> Label<N>& = default;
    auto operator=(Label<N>&& src) noexcept -> Label<N>& = default;
    ~Label() = default;

    void set_f(double _f) { f = _f; }

    void set_job(Job* _job) { label_job = _job; };

    void set_node_id(NodeId _id) { node_id = _id; }

    void reset() {
        f = std::numeric_limits<double>::max();
        prev_label = nullptr;
        high = false;
    }

    [[nodiscard]] auto get_f() const -> double { return f; }

    auto get_f() -> double& { return f; }

    auto get_previous() -> Label<N>* { return prev_label; }

    auto prev_job_backward() -> Job* { return prev_job; }

    auto get_high() -> bool { return high; }

    auto get_job() -> Job* { return label_job; }

    auto get_weight() -> int { return weight; }

    auto get_node_id() -> NodeId& { return node_id; }

    auto prev_job_forward() -> Job* {
        return get_previous() == nullptr ? nullptr : get_previous()->get_job();
    }

    void forward_update(double _f, Label<N>* _prev, bool _high = false) {
        f = _f;
        prev_label = _prev;
        high = _high;
    }

    void forward_update(double _f, Label<N>& _node) {
        f = _f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void forward_update(Label<N>& _node) {
        f = _node.f;
        prev_label = _node.prev_label;
        high = _node.high;
    }

    void backward_update(double _f, bool _high = false) {
        f = _f;
        high = _high;
    }

    void backward_update(Label<N>* _n, double _f = .0, bool _high = false) {
        f = _f;
        prev_job = _high ? get_job() : _n->prev_job;
        high = _high;
        prev_label = _n;
    }
};

#endif  // LABEL_HPP
