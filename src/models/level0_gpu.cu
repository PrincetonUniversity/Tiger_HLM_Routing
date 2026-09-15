#include "level0_gpu.hpp"

#include <thrust/device_vector.h>
#include <thrust/for_each.h>
#include <thrust/copy.h>
#include <thrust/iterator/counting_iterator.h>

#include <cuda_runtime.h>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

// =============================================================================
// Fused RK4: updates state AND stores to batch output buffer on GPU
// =============================================================================
struct FusedRK4BatchFunctor {
    double*       __restrict__ q;
    float*        __restrict__ batch_out;
    const double* __restrict__ A_h;
    const double* __restrict__ lambda_1;
    const double* __restrict__ invtau;
    const float*  __restrict__ runoff;
    size_t n0, nTime, resolution, batch_size, step_in_batch;
    double dt, t;

    __device__ double rhs(double qv, double r, double ah, double lam, double itau) const {
        double qs = qv < 1e-8 ? 1e-8 : qv;
        double rt = (r * 0.001 / 60.0) * ah / 60.0;
        return itau * pow(qs, lam) * (-qs + rt);
    }

    __device__ void operator()(size_t k) const {
        double qk = q[k];

        // Evaluate runoff at each RK4 intermediate time, matching Boost's CPU behavior
        auto get_runoff = [&](double t_eval) -> double {
            size_t ri = static_cast<size_t>(t_eval) / resolution;
            if (ri >= nTime) ri = nTime - 1;
            return static_cast<double>(runoff[ri * n0 + k]);
        };

        double ah   = A_h[k];
        double lam  = lambda_1[k];
        double itau = invtau[k];

        double r1 = get_runoff(t);            // t
        double r2 = get_runoff(t + 0.5*dt);   // t + dt/2 (for k2 and k3)
        double r4 = get_runoff(t + dt);        // t + dt   (for k4)

        double k1 = rhs(qk,              r1, ah, lam, itau);
        double k2 = rhs(qk + 0.5*dt*k1, r2, ah, lam, itau);
        double k3 = rhs(qk + 0.5*dt*k2, r2, ah, lam, itau);
        double k4 = rhs(qk + dt*k3,     r4, ah, lam, itau);

        double q_new = qk + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        if (q_new < 1e-8) q_new = 1e-8;
        q[k] = q_new;
        batch_out[k * batch_size + step_in_batch] = static_cast<float>(q_new);
    }
};

// =============================================================================
// Persistent context
// =============================================================================
namespace {

struct Level0Ctx {
    bool   built = false;
    size_t n0    = 0;
    std::vector<size_t> node_index_h;
    std::vector<size_t> local_index_h;
    std::vector<int>    stream_id_h;
    thrust::device_vector<double> A_h_d, lambda_1_d, invtau_d;
};

Level0Ctx g_ctx;

void ensureCtx(const ModelSetup& setup,
               const Partition& part,
               const std::vector<size_t>& nodes)
{
    std::vector<size_t> owned;
    for (size_t idx : nodes) {
        if (part.owns(idx)) owned.push_back(idx);
    }
    if (g_ctx.built && g_ctx.n0 == owned.size()) return;
    size_t n0 = owned.size();
    std::vector<size_t> gidx(n0), lidx(n0);
    std::vector<int>    sid(n0);
    std::vector<double> ah(n0), lam(n0), itau(n0);
    for (size_t k = 0; k < n0; ++k) {
        const NodeInfo& nd = setup.node_map.at(owned[k]);
        gidx[k] = nd.index;
        lidx[k] = part.local_of[nd.index];
        sid[k]  = nd.stream_id;
        ah[k]   = nd.params[0];
        lam[k]  = nd.params[2];
        itau[k] = 60.0 * nd.params[3] / ((1.0 - lam[k]) * nd.params[1]);
    }
    g_ctx.built         = true;
    g_ctx.n0            = n0;
    g_ctx.node_index_h  = gidx;
    g_ctx.local_index_h = lidx;
    g_ctx.stream_id_h   = sid;
    g_ctx.A_h_d         = ah;
    g_ctx.lambda_1_d    = lam;
    g_ctx.invtau_d      = itau;
}

size_t pickBatchSize(size_t n0, size_t n_steps) {
    const size_t GPU_BUDGET = 74ULL * 1024 * 1024 * 1024;
    size_t per_step = n0 * sizeof(float);
    size_t max_batch = GPU_BUDGET / per_step;
    if (max_batch < 1) max_batch = 1;
    if (max_batch > n_steps) max_batch = n_steps;
    return max_batch;
}

} // anonymous namespace

// =============================================================================
// Public entry point — batched GPU integration
// =============================================================================
void IntegrateLevel0GPU(const ModelSetup& setup,
                        const Partition& part,
                        const RunoffData& runoff,
                        std::vector<float>& results,
                        const std::vector<size_t>& nodes_at_level,
                        size_t n_steps,
                        size_t tc,
                        std::vector<float>& q_final)
{
    ensureCtx(setup, part, nodes_at_level);
    const size_t n0 = g_ctx.n0;
    if (n0 == 0) return;

    auto wall_start = std::chrono::high_resolution_clock::now();

    // ======== GATHER ========
    std::vector<float> h_ro(runoff.nTime * n0);
    #pragma omp parallel for
    for (size_t k = 0; k < n0; ++k) {
        auto it = runoff.idToIndex.find(g_ctx.stream_id_h[k]);
        if (it == runoff.idToIndex.end()) {
            fprintf(stderr, "IntegrateLevel0GPU: stream_id %d missing\n",
                    g_ctx.stream_id_h[k]);
            exit(EXIT_FAILURE);
        }
        const float* src = &runoff.data[it->second * runoff.nTime];
        for (size_t t = 0; t < runoff.nTime; ++t)
            h_ro[t * n0 + k] = src[t];
    }

    std::vector<double> h_q0(n0);
    if (tc > 0) {
        std::cout << "    [GPU-L0 READ DEBUG] tc=" << tc << " reading q_final for chunk start:\n";
    }
    for (size_t k = 0; k < n0; ++k) {
        if (tc == 0) {
            h_q0[k] = setup.uini(g_ctx.stream_id_h[k]);
        } else {
            size_t local_idx = g_ctx.local_index_h[k];
            h_q0[k] = q_final[local_idx];
            if (k < 5) {
                std::cout << "      k=" << k
                          << " global=" << g_ctx.node_index_h[k]
                          << " local=" << local_idx
                          << " stream_id=" << g_ctx.stream_id_h[k]
                          << " q_final_read=" << h_q0[k]
                          << "\n";
            }
            if (h_q0[k] <= 0.0) {
                fprintf(stderr, "Level0GPU: non-positive q0 at link %zu\n",
                        g_ctx.node_index_h[k]);
                exit(EXIT_FAILURE);
            }
        }
    }

    thrust::device_vector<double> d_q      = h_q0;
    thrust::device_vector<float>  d_runoff = h_ro;

    size_t batch_size = pickBatchSize(n0, n_steps - 1);
    thrust::device_vector<float> d_batch(n0 * batch_size);
    std::vector<float> h_batch(n0 * batch_size);

    cudaDeviceSynchronize();
    auto gather_end = std::chrono::high_resolution_clock::now();

    for (size_t k = 0; k < n0; ++k) {
        size_t local = g_ctx.local_index_h[k];
        results[local * n_steps + 0] = static_cast<float>(std::max(h_q0[k], 1e-8));
    }

    // ======== INTEGRATE IN BATCHES ========
    auto integrate_start = std::chrono::high_resolution_clock::now();

    double* q_ptr    = thrust::raw_pointer_cast(d_q.data());
    double* ah_ptr   = thrust::raw_pointer_cast(g_ctx.A_h_d.data());
    double* lam_ptr  = thrust::raw_pointer_cast(g_ctx.lambda_1_d.data());
    double* itau_ptr = thrust::raw_pointer_cast(g_ctx.invtau_d.data());
    float*  ro_ptr   = thrust::raw_pointer_cast(d_runoff.data());
    float*  bout_ptr = thrust::raw_pointer_cast(d_batch.data());

    size_t resolution = static_cast<size_t>(setup.config.runoff_resolution);
    double dt = setup.config.dt;

    size_t total_steps = n_steps - 1;
    size_t steps_done = 0;
    size_t n_batches = 0;
    double copy_ms = 0, scatter_ms = 0;

    while (steps_done < total_steps) {
        size_t this_batch = std::min(batch_size, total_steps - steps_done);

        for (size_t b = 0; b < this_batch; ++b) {
            size_t global_step = steps_done + b;
            double t = global_step * dt;

            thrust::for_each(
                thrust::counting_iterator<size_t>(0),
                thrust::counting_iterator<size_t>(n0),
                FusedRK4BatchFunctor{
                    q_ptr, bout_ptr,
                    ah_ptr, lam_ptr, itau_ptr, ro_ptr,
                    n0, runoff.nTime, resolution, this_batch, b,
                    dt, t
                });
        }

        cudaDeviceSynchronize();
        auto copy_start = std::chrono::high_resolution_clock::now();

        cudaMemcpy(h_batch.data(), bout_ptr,
                   n0 * this_batch * sizeof(float),
                   cudaMemcpyDeviceToHost);

        auto copy_end = std::chrono::high_resolution_clock::now();
        copy_ms += std::chrono::duration<double, std::milli>(copy_end - copy_start).count();

        auto scatter_start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for
        for (size_t k = 0; k < n0; ++k) {
            size_t local = g_ctx.local_index_h[k];
            for (size_t b = 0; b < this_batch; ++b) {
                size_t out_step = steps_done + b + 1;
                results[local * n_steps + out_step] = h_batch[k * this_batch + b];
            }
        }

        auto scatter_end = std::chrono::high_resolution_clock::now();
        scatter_ms += std::chrono::duration<double, std::milli>(scatter_end - scatter_start).count();

        steps_done += this_batch;
        n_batches++;
    }
    // ======== DEBUG: dump GPU state for chunk handoff ========
    {
        std::vector<double> h_q_final_gpu(n0);
        cudaMemcpy(h_q_final_gpu.data(), q_ptr, n0 * sizeof(double), cudaMemcpyDeviceToHost);

        size_t last_step = n_steps - 1;
        std::cout << "    [GPU-L0 DEBUG] tc=" << tc << " n0=" << n0
                  << " n_steps=" << n_steps << " last_step=" << last_step << "\n";

        // Show first 5 links: device q vs results buffer vs q_final handoff
        size_t show = std::min(n0, (size_t)5);
        for (size_t k = 0; k < show; ++k) {
            size_t local = g_ctx.local_index_h[k];
            size_t gidx  = g_ctx.node_index_h[k];
            int    sid   = g_ctx.stream_id_h[k];
            float  res_first = results[local * n_steps + 0];
            float  res_last  = results[local * n_steps + last_step];
            std::cout << "      k=" << k
                      << " global=" << gidx
                      << " local=" << local
                      << " stream_id=" << sid
                      << " q0_input=" << h_q0[k]
                      << " d_q_final=" << h_q_final_gpu[k]
                      << " results[0]=" << res_first
                      << " results[last]=" << res_last
                      << "\n";
        }

        // Show last 3 links too
        for (size_t k = (n0 > 5 ? n0 - 3 : 5); k < n0; ++k) {
            size_t local = g_ctx.local_index_h[k];
            float  res_last  = results[local * n_steps + last_step];
            std::cout << "      k=" << k
                      << " global=" << g_ctx.node_index_h[k]
                      << " local=" << local
                      << " stream_id=" << g_ctx.stream_id_h[k]
                      << " q0_input=" << h_q0[k]
                      << " d_q_final=" << h_q_final_gpu[k]
                      << " results[last]=" << res_last
                      << "\n";
        }
        std::cout << std::flush;
    }

    cudaDeviceSynchronize();
    auto integrate_end = std::chrono::high_resolution_clock::now();

    // ======== TIMING ========
    double gather_ms    = std::chrono::duration<double, std::milli>(gather_end      - wall_start).count();
    double integrate_ms = std::chrono::duration<double, std::milli>(integrate_end   - integrate_start).count();
    double total_ms     = std::chrono::duration<double, std::milli>(integrate_end   - wall_start).count();
    double kernel_ms    = integrate_ms - copy_ms - scatter_ms;

    double gpu_gb = (n0 * sizeof(double) +
                     runoff.nTime * n0 * sizeof(float) +
                     n0 * 3 * sizeof(double) +
                     n0 * batch_size * sizeof(float)) / 1e9;

    std::cout << "    [GPU-L0 batched] links=" << n0
              << " steps=" << n_steps
              << " batch=" << batch_size
              << " batches=" << n_batches
              << " gpu_mem=" << gpu_gb << "GB"
              << "\n        gather=" << gather_ms << "ms"
              << "  kernels=" << kernel_ms << "ms"
              << "  copy=" << copy_ms << "ms"
              << "  scatter=" << scatter_ms << "ms"
              << "  total=" << total_ms << "ms"
              << std::endl;
}
void FreeLevel0GPU() {
    g_ctx.A_h_d.clear();
    g_ctx.A_h_d.shrink_to_fit();
    g_ctx.lambda_1_d.clear();
    g_ctx.lambda_1_d.shrink_to_fit();
    g_ctx.invtau_d.clear();
    g_ctx.invtau_d.shrink_to_fit();
    g_ctx.node_index_h.clear();
    g_ctx.local_index_h.clear();
    g_ctx.stream_id_h.clear();
    g_ctx.built = false;
    g_ctx.n0    = 0;
}
