#include "level0_gpu.hpp"

// Thrust
#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>

// Boost.Odeint with its built-in Thrust backend
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>

#include <cuda_runtime.h>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace boost::numeric::odeint;

// =============================================================================
// State type + stepper: use exactly what boost provides for Thrust
// =============================================================================
typedef thrust::device_vector<double> state_type;

typedef runge_kutta4
    state_type,     // state
    double,         // value
    state_type,     // deriv
    double,         // time
    thrust_algebra,
    thrust_operations
> rk4_thrust_type;

// =============================================================================
// Persistent context — rebuilt only if level-0 set changes
// =============================================================================
namespace {

struct Level0Ctx {
    bool   built = false;
    size_t n0    = 0;
    std::vector<size_t> node_index_h;
    std::vector<int>    stream_id_h;
    thrust::device_vector<double> A_h_d, lambda_1_d, invtau_d;
};

Level0Ctx g_ctx;

void ensureCtx(const ModelSetup& setup,
               const std::vector<size_t>& nodes)
{
    if (g_ctx.built && g_ctx.n0 == nodes.size()) return;

    size_t n0 = nodes.size();
    std::vector<size_t> idx(n0);
    std::vector<int>    sid(n0);
    std::vector<double> ah(n0), lam(n0), itau(n0);

    for (size_t k = 0; k < n0; ++k) {
        const NodeInfo& nd = setup.node_map.at(nodes[k]);
        idx[k]  = nd.index;
        sid[k]  = nd.stream_id;
        ah[k]   = nd.params[0];
        lam[k]  = nd.params[2];
        itau[k] = 60.0 * nd.params[3] / ((1.0 - lam[k]) * nd.params[1]);
    }

    g_ctx.built        = true;
    g_ctx.n0           = n0;
    g_ctx.node_index_h = idx;
    g_ctx.stream_id_h  = sid;
    g_ctx.A_h_d        = ah;
    g_ctx.lambda_1_d   = lam;
    g_ctx.invtau_d     = itau;
}

// =============================================================================
// Device functor: computes dQ/dt for one link
// Called from OdeSystem via thrust::for_each — this is the ONLY device code
// =============================================================================
struct RHSKernel {
    const double* __restrict__ q;
    double*       __restrict__ dqdt;
    const double* __restrict__ A_h;
    const double* __restrict__ lambda_1;
    const double* __restrict__ invtau;
    const float*  __restrict__ runoff;   // time-major: [t * n0 + k]
    size_t n0, nTime, resolution;
    double t;

    __device__ void operator()(size_t k) const {
        double qv = q[k];
        double qs = qv < 1e-8 ? 1e-8 : qv;
        size_t ri = static_cast<size_t>(t) / resolution;
        if (ri >= nTime) ri = nTime - 1;
        double r  = static_cast<double>(runoff[ri * n0 + k]);
        double rt = (r * 0.001 / 60.0) * A_h[k] / 60.0;
        dqdt[k]   = invtau[k] * pow(qs, lambda_1[k]) * (-qv + rt);
    }
};

// =============================================================================
// System callable — this is what odeint calls at each RK4 stage
// Host-side function that dispatches RHSKernel via thrust::for_each
// Each call = 1 kernel launch
// =============================================================================
struct OdeSystem {
    const double* A_h;
    const double* lambda_1;
    const double* invtau;
    const float*  runoff;
    size_t n0, nTime, resolution;

    void operator()(const state_type& x, state_type& dxdt, double t) const {
        thrust::for_each(
            thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(n0),
            RHSKernel{
                thrust::raw_pointer_cast(x.data()),
                thrust::raw_pointer_cast(dxdt.data()),
                A_h, lambda_1, invtau, runoff,
                n0, nTime, resolution, t
            });
    }
};

// =============================================================================
// Observer — odeint calls this at each output step (host-side)
// Stores device state into a device buffer via thrust::for_each
// Each call = 1 kernel launch
// =============================================================================
struct StoreKernel {
    const double* __restrict__ q;
    float*        __restrict__ out;
    size_t n_steps, step;

    __device__ void operator()(size_t k) const {
        out[k * n_steps + step] = fmaxf(static_cast<float>(q[k]), 1e-8f);
    }
};

struct Observer {
    float* d_out;
    size_t n0, n_steps;
    double dt;

    void operator()(const state_type& x, double t) {
        size_t s = static_cast<size_t>(t / dt);
        if (s >= n_steps) s = n_steps - 1;
        thrust::for_each(
            thrust::device,
            thrust::counting_iterator<size_t>(0),
            thrust::counting_iterator<size_t>(n0),
            StoreKernel{
                thrust::raw_pointer_cast(x.data()),
                d_out, n_steps, s
            });
    }
};

} // anonymous namespace

// =============================================================================
// Public entry point
// =============================================================================
void IntegrateLevel0GPU(const ModelSetup& setup,
                        const RunoffData& runoff,
                        std::vector<float>& results,
                        const std::vector<size_t>& nodes_at_level,
                        size_t n_steps,
                        size_t tc,
                        std::vector<float>& q_final)
{
    ensureCtx(setup, nodes_at_level);
    const size_t n0 = g_ctx.n0;
    if (n0 == 0) return;

    auto wall_start = std::chrono::high_resolution_clock::now();

    // ======== GATHER: transpose runoff, gather q0, upload to device ========
    auto gather_start = std::chrono::high_resolution_clock::now();

    // Transpose runoff from link-major [link * nTime + t] to time-major [t * n0 + k]
    // so that at each RK4 stage, consecutive threads read consecutive memory
    std::vector<float> h_ro(runoff.nTime * n0);
    #pragma omp parallel for
    for (size_t k = 0; k < n0; ++k) {
        auto it = runoff.idToIndex.find(g_ctx.stream_id_h[k]);
        if (it == runoff.idToIndex.end()) {
            fprintf(stderr, "IntegrateLevel0GPU: stream_id %d missing from runoff\n",
                    g_ctx.stream_id_h[k]);
            exit(EXIT_FAILURE);
        }
        const float* src = &runoff.data[it->second * runoff.nTime];
        for (size_t t = 0; t < runoff.nTime; ++t)
            h_ro[t * n0 + k] = src[t];
    }
    thrust::device_vector<float> d_runoff = h_ro;

    // Initial conditions — same tc==0 / tc>0 branch as CPU routing.cpp
    std::vector<double> h_q0(n0);
    for (size_t k = 0; k < n0; ++k) {
        if (tc == 0) {
            h_q0[k] = setup.uini(g_ctx.stream_id_h[k]);
        } else {
            h_q0[k] = q_final[g_ctx.node_index_h[k]];
            if (h_q0[k] <= 0.0) {
                fprintf(stderr, "Level0GPU: non-positive q0 at link %zu\n",
                        g_ctx.node_index_h[k]);
                exit(EXIT_FAILURE);
            }
        }
    }
    state_type d_state = h_q0;  // host→device upload

    // Device output buffer: [n0 * n_steps], link-major (same layout as CPU results)
    thrust::device_vector<float> d_out(n0 * n_steps);

    cudaDeviceSynchronize();
    auto gather_end = std::chrono::high_resolution_clock::now();

    // ======== INTEGRATE: boost odeint RK4 with its built-in Thrust backend ========
    auto integrate_start = std::chrono::high_resolution_clock::now();

    OdeSystem sys{
        thrust::raw_pointer_cast(g_ctx.A_h_d.data()),
        thrust::raw_pointer_cast(g_ctx.lambda_1_d.data()),
        thrust::raw_pointer_cast(g_ctx.invtau_d.data()),
        thrust::raw_pointer_cast(d_runoff.data()),
        n0, runoff.nTime,
        static_cast<size_t>(setup.config.runoff_resolution)
    };

    Observer obs{
        thrust::raw_pointer_cast(d_out.data()),
        n0, n_steps, setup.config.dt
    };

    // This is the same integrate_const you already use on CPU,
    // just with a thrust-backed stepper instead of a scalar one.
    // Internally: each RK4 step = 4 system calls + ~4 algebra ops + 1 observer
    //           = ~9 kernel launches per step
    rk4_thrust_type stepper;
    double t0 = 0.0;
    double t1 = (n_steps - 1) * setup.config.dt;
    integrate_const(stepper, sys, d_state, t0, t1, setup.config.dt, obs);

    cudaDeviceSynchronize();
    auto integrate_end = std::chrono::high_resolution_clock::now();

    // ======== SCATTER: copy device results back into host `results` buffer ========
    auto scatter_start = std::chrono::high_resolution_clock::now();

    std::vector<float> h_out(n0 * n_steps);
    thrust::copy(d_out.begin(), d_out.end(), h_out.begin());

    #pragma omp parallel for
    for (size_t k = 0; k < n0; ++k) {
        size_t row = g_ctx.node_index_h[k];
        std::copy(h_out.data() + k * n_steps,
                  h_out.data() + (k + 1) * n_steps,
                  results.data() + row * n_steps);
    }

    auto scatter_end = std::chrono::high_resolution_clock::now();
    auto wall_end    = std::chrono::high_resolution_clock::now();

    // ======== TIMING REPORT ========
    double gather_ms    = std::chrono::duration<double, std::milli>(gather_end    - gather_start).count();
    double integrate_ms = std::chrono::duration<double, std::milli>(integrate_end - integrate_start).count();
    double scatter_ms   = std::chrono::duration<double, std::milli>(scatter_end   - scatter_start).count();
    double total_ms     = std::chrono::duration<double, std::milli>(wall_end      - wall_start).count();

    // ~9 kernel launches per step: 4 system evals + ~4 thrust_algebra ops + 1 observer
    size_t est_launches = 9 * (n_steps > 0 ? n_steps - 1 : 0) + 1;

    std::cout << "    [GPU-L0 odeint+thrust] links=" << n0
              << " steps=" << n_steps
              << "\n        gather="    << gather_ms    << "ms"
              << "  integrate=" << integrate_ms << "ms"
              << "  scatter="   << scatter_ms   << "ms"
              << "  total="     << total_ms     << "ms"
              << "\n        est_kernel_launches=" << est_launches
              << "  (overhead_per_launch~="
              << (est_launches > 0 ? integrate_ms / est_launches : 0.0) << "ms)"
              << std::endl;
}