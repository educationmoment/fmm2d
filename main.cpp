#include "fmm.hpp"
#include <chrono>
#include <cstdio>
#include <fstream>

using namespace fmm;

struct Err { double linf, l1; };


static Err pointSource(int n, bool exactInit) {
    const double h = 2.0 / (n - 1);
    FastMarching2D s(n, n, h, h);
    const int c = n / 2;
    const double r0 = 0.1;
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const double x = -1 + i * h, y = -1 + j * h, r = std::hypot(x, y);
            if ((exactInit && r <= r0) || (i == c && j == c)) s.addSource(i, j, r);
        }
    s.solve();
    Err e{0, 0};
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const double x = -1 + i * h, y = -1 + j * h;
            const double d = std::abs(s.at(i, j) - std::hypot(x, y));
            e.linf = std::max(e.linf, d);
            e.l1 += d * h * h;
        }
    return e;
}

static void study(bool exactInit) {
    std::printf("\n%s\n", exactInit ? "Point source, exact init in r<=0.1 disc"
                                      : "Point source, single-node init");
    std::printf("%6s %12s %7s %12s %7s\n", "N", "Linf", "rate", "L1", "rate");
    Err prev{};
    for (int k = 0, n = 51; k < 5; ++k, n = 2 * n - 1) {
        Err e = pointSource(n, exactInit);
        if (k == 0) std::printf("%6d %12.4e %7s %12.4e %7s\n", n, e.linf, "-", e.l1, "-");
        else std::printf("%6d %12.4e %7.3f %12.4e %7.3f\n", n, e.linf,
                         std::log2(prev.linf / e.linf), e.l1, std::log2(prev.l1 / e.l1));
        prev = e;
    }
}

int main() {
    study(false);
    study(true);

    // obstaclez wall at x = 0.5 with a slit source at bottom-left.
    const int n = 201;
    const double h = 1.0 / (n - 1);
    std::vector<double> F(n * n, 1.0);
    for (int j = 0; j < n; ++j)
        if (std::abs(j * h - 0.5) > 0.05) F[j * n + n / 2] = 0.0;
    FastMarching2D s(n, n, h, h);
    s.setSpeed(F);
    s.addSource(0, 0);
    auto t0 = std::chrono::steady_clock::now();
    s.solve();
    auto t1 = std::chrono::steady_clock::now();
    std::printf("\nObstacle test %dx%d: T(top-right)=%.4f (straight line would be %.4f), %.2f ms\n",
                n, n, s.at(n - 1, n - 1), std::sqrt(2.0),
                std::chrono::duration<double, std::milli>(t1 - t0).count());
    std::ofstream out("obstacle_T.csv");
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) out << s.at(i, j) << (i + 1 < n ? "," : "\n");
    }

    const int N = 2001;
    FastMarching2D big(N, N, 1.0 / (N - 1), 1.0 / (N - 1));
    big.addSource(N / 2, N / 2);
    t0 = std::chrono::steady_clock::now();
    big.solve();
    t1 = std::chrono::steady_clock::now();
    std::printf("Timing %dx%d (%.1f M nodes): %.1f ms\n", N, N, N * N / 1e6,
                std::chrono::duration<double, std::milli>(t1 - t0).count());
}
