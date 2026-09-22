#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace fmm {

inline constexpr double kInf = std::numeric_limits<double>::infinity();

// Sethian Sec. 8.2: O(log N)
class IndexedMinHeap {
public:
    IndexedMinHeap(const std::vector<double>& key, std::size_t n)
        : key_(key), pos_(n, -1) { heap_.reserve(n / 8 + 16); }

    bool empty() const { return heap_.empty(); }
    bool contains(std::size_t i) const { return pos_[i] >= 0; }

    void push(std::size_t i) {
        heap_.push_back(i);
        pos_[i] = static_cast<std::ptrdiff_t>(heap_.size() - 1);
        siftUp(heap_.size() - 1);
    }
    void decrease(std::size_t i) { siftUp(static_cast<std::size_t>(pos_[i])); }

    std::size_t pop() {
        const std::size_t top = heap_.front();
        swapNodes(0, heap_.size() - 1);
        heap_.pop_back();
        pos_[top] = -1;
        if (!heap_.empty()) siftDown(0);
        return top;
    }

private:
    const std::vector<double>& key_;
    std::vector<std::size_t> heap_;
    std::vector<std::ptrdiff_t> pos_;

    bool less(std::size_t a, std::size_t b) const { return key_[heap_[a]] < key_[heap_[b]]; }
    void swapNodes(std::size_t a, std::size_t b) {
        std::swap(heap_[a], heap_[b]);
        pos_[heap_[a]] = static_cast<std::ptrdiff_t>(a);
        pos_[heap_[b]] = static_cast<std::ptrdiff_t>(b);
    }
    void siftUp(std::size_t k) {
        while (k > 0) {
            const std::size_t p = (k - 1) / 2;
            if (!less(k, p)) break;
            swapNodes(k, p);
            k = p;
        }
    }
    void siftDown(std::size_t k) {
        const std::size_t n = heap_.size();
        for (;;) {
            const std::size_t l = 2 * k + 1;
            if (l >= n) break;
            std::size_t m = l;
            if (l + 1 < n && less(l + 1, l)) m = l + 1;
            if (!less(m, k)) break;
            swapNodes(k, m);
            k = m;
        }
    }
};

enum class State : std::uint8_t { Far, Trial, Known };

class FastMarching2D {
public:
    FastMarching2D(int nx, int ny, double dx, double dy)
        : nx_(nx), ny_(ny), dx_(dx), dy_(dy),
          F_(static_cast<std::size_t>(nx) * ny, 1.0),
          T_(F_.size(), kInf), state_(F_.size(), State::Far) {
        if (nx < 1 || ny < 1 || dx <= 0 || dy <= 0)
            throw std::invalid_argument("bad grid");
    }

    // F <= 0 marks an obstacle: T stays +inf there.
    void setSpeed(std::vector<double> F) {
        if (F.size() != F_.size()) throw std::invalid_argument("speed size mismatch");
        F_ = std::move(F);
    }

    void addSource(int i, int j, double t0 = 0.0) {
        if (!inside(i, j)) throw std::out_of_range("source outside grid");
        const std::size_t k = idx(i, j);
        T_[k] = std::min(T_[k], t0);
        state_[k] = State::Known;
    }

    void solve() {
        IndexedMinHeap heap(T_, T_.size());   // T_ must not reallocate from here on

        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
                if (state_[idx(i, j)] == State::Known) relaxNeighbours(i, j, heap);

        while (!heap.empty()) {
            const std::size_t k = heap.pop();
            state_[k] = State::Known;
            relaxNeighbours(static_cast<int>(k % nx_), static_cast<int>(k / nx_), heap);
        }
    }

    double at(int i, int j) const { return T_[idx(i, j)]; }
    const std::vector<double>& arrival() const { return T_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }

private:
    int nx_, ny_;
    double dx_, dy_;
    std::vector<double> F_, T_;
    std::vector<State> state_;

    std::size_t idx(int i, int j) const { return static_cast<std::size_t>(j) * nx_ + i; }
    bool inside(int i, int j) const { return i >= 0 && i < nx_ && j >= 0 && j < ny_; }

    double knownValue(int i, int j) const {
        if (!inside(i, j)) return kInf;
        const std::size_t k = idx(i, j);
        return state_[k] == State::Known ? T_[k] : kInf;
    }

    double localSolve(int i, int j) const {
        const double f = F_[idx(i, j)];
        double t0 = std::min(knownValue(i - 1, j), knownValue(i + 1, j));
        double t1 = std::min(knownValue(i, j - 1), knownValue(i, j + 1));
        double h0 = dx_, h1 = dy_;
        if (t1 < t0) { std::swap(t0, t1); std::swap(h0, h1); }
        if (t0 == kInf) return kInf;

        double T = t0 + h0 / f;

        if (T > t1) {
            const double w0 = 1.0 / (h0 * h0), w1 = 1.0 / (h1 * h1);
            const double d = t0 - t1;
            const double disc = (w0 + w1) / (f * f) - w0 * w1 * d * d;
            if (disc >= 0.0) T = (w0 * t0 + w1 * t1 + std::sqrt(disc)) / (w0 + w1);
        }
        return T;
    }

    void relaxNeighbours(int i, int j, IndexedMinHeap& heap) {
        static constexpr int di[4] = {-1, 1, 0, 0};
        static constexpr int dj[4] = {0, 0, -1, 1};
        for (int n = 0; n < 4; ++n) {
            const int a = i + di[n], b = j + dj[n];
            if (!inside(a, b)) continue;
            const std::size_t k = idx(a, b);
            if (state_[k] == State::Known || F_[k] <= 0.0) continue;
            const double t = localSolve(a, b);
            if (t < T_[k]) {
                T_[k] = t;
                if (state_[k] == State::Trial) heap.decrease(k);
                else { state_[k] = State::Trial; heap.push(k); }
            }
        }
    }
};

}  // namespace fmm
