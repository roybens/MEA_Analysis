/**
 * mxw_streamer.cpp
 * ================
 * Taps the MaxOne raw data stream from mxwserver and writes selected
 * channels into a POSIX shared memory ring buffer for the Python
 * oscilloscope to read.
 *
 * This runs ALONGSIDE MaxLab Live Scope — it does not interfere with
 * Scope or with setup.py. It only reads the data stream.
 *
 * Compile
 * -------
 *   g++ -std=c++17 -O3 -o mxw_streamer mxw_streamer.cpp \
 *       -I/path/to/maxlab_lib/include \
 *       -L/path/to/maxlab_lib/lib \
 *       -lmaxlab -lrt -lpthread \
 *       -Wl,-rpath,/path/to/maxlab_lib/lib
 *
 * Run
 * ---
 *   ./mxw_streamer <ch1> <ch2> <ch3> [ch4]
 *
 *   Channel numbers come from MaxLab Live Scope GUI — hover over an
 *   electrode and Scope shows its channel number (0–1023).
 *
 * Example
 * -------
 *   ./mxw_streamer 42 55 67 88
 */

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <atomic>
#include <thread>
#include <chrono>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "maxlab/maxlab.h"

// ─────────────────────────────────────────────────────────────────────────────
// Shared memory constants — must match oscilloscope.py exactly
// ─────────────────────────────────────────────────────────────────────────────

static constexpr const char* SHM_NAME     = "/mxw_scope";
static constexpr uint32_t    SHM_MAGIC    = 0xC0FFEE42;
static constexpr int         MAX_CHANNELS = 8;
static constexpr int         RING_FRAMES  = 40000;   // 2 s @ 20 kHz
static constexpr size_t      DATA_OFFSET  = 4096;    // page-aligned start of data
static constexpr size_t      DATA_SIZE    = RING_FRAMES * MAX_CHANNELS * sizeof(float);
static constexpr size_t      SHM_TOTAL    = DATA_OFFSET + DATA_SIZE;

static constexpr uint8_t     TARGET_WELL  = 0;       // MaxOne = 0

// ─────────────────────────────────────────────────────────────────────────────
// Shared memory header — Python reads this to know layout and current write pos
// Fields must be in this exact order (matches HEADER_FMT in oscilloscope.py)
// ─────────────────────────────────────────────────────────────────────────────

struct ShmHeader {
    uint32_t magic;                      // set last — Python polls this to know data is ready
    uint32_t n_channels;                 // number of active channels
    uint16_t channel_ids[MAX_CHANNELS];  // which channel indices are stored
    uint64_t write_pos;                  // absolute frame counter (not modded)
    uint64_t frame_number;               // hardware frame number of latest write
    uint32_t sample_rate;                // Hz
    float    lsb_uv;                     // µV per raw ADC unit
};

// ─────────────────────────────────────────────────────────────────────────────

static std::atomic<bool> g_running{true};
void on_signal(int) { g_running = false; }

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <ch1> [ch2] [ch3] [ch4]\n"
                  << "  Channel numbers from MaxLab Live Scope GUI (0–1023)\n"
                  << "  Max " << MAX_CHANNELS << " channels\n";
        return 1;
    }

    // ── Parse channel indices ─────────────────────────────────────────────
    std::vector<int> channels;
    for (int i = 1; i < argc && (int)channels.size() < MAX_CHANNELS; ++i)
        channels.push_back(std::atoi(argv[i]));
    const int n_ch = (int)channels.size();

    std::cout << "[streamer] Channels to stream: ";
    for (int c : channels) std::cout << c << " ";
    std::cout << "\n";

    // ── Create shared memory ──────────────────────────────────────────────
    shm_unlink(SHM_NAME);  // remove any stale segment from a previous run

    int shm_fd = shm_open(SHM_NAME, O_CREAT | O_RDWR, 0666);
    if (shm_fd < 0) { perror("shm_open failed"); return 1; }

    if (ftruncate(shm_fd, SHM_TOTAL) < 0) { perror("ftruncate failed"); return 1; }

    void* shm_ptr = mmap(nullptr, SHM_TOTAL,
                          PROT_READ | PROT_WRITE, MAP_SHARED,
                          shm_fd, 0);
    if (shm_ptr == MAP_FAILED) { perror("mmap failed"); return 1; }
    close(shm_fd);

    memset(shm_ptr, 0, SHM_TOTAL);

    // ── Fill header (magic written last so Python knows it's ready) ───────
    ShmHeader* hdr     = reinterpret_cast<ShmHeader*>(shm_ptr);
    float*     data    = reinterpret_cast<float*>((char*)shm_ptr + DATA_OFFSET);

    hdr->n_channels  = (uint32_t)n_ch;
    hdr->sample_rate = 20000;
    hdr->lsb_uv      = 2.34f;   // µV/LSB at gain=512 (approximate)
    hdr->write_pos   = 0;
    hdr->frame_number = 0;
    for (int i = 0; i < n_ch; ++i)
        hdr->channel_ids[i] = (uint16_t)channels[i];

    __sync_synchronize();   // ensure all header writes land before magic
    hdr->magic = SHM_MAGIC; // Python starts reading only after seeing this

    std::cout << "[streamer] Shared memory ready: " << SHM_NAME
              << "  total=" << (SHM_TOTAL / 1024) << " KB\n";

    // ── Open raw stream ───────────────────────────────────────────────────
    maxlab::checkVersions();
    maxlab::verifyStatus(maxlab::DataStreamerRaw_open());
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    std::cout << "[streamer] Raw stream open. Streaming to shared memory...\n";

    signal(SIGINT,  on_signal);
    signal(SIGTERM, on_signal);

    maxlab::RawFrameData frame;
    uint64_t frames_written = 0;

    // ── Main loop — as fast as frames arrive (~20 kHz) ────────────────────
    while (g_running) {
        maxlab::Status st = maxlab::DataStreamerRaw_receiveNextFrame(&frame);
        if (st == maxlab::Status::MAXLAB_NO_FRAME) continue;
        if (frame.frameInfo.well_id != TARGET_WELL) continue;
        if (frame.frameInfo.corrupted) continue;

        // Write into ring buffer slot
        uint64_t slot = hdr->write_pos % RING_FRAMES;
        float*   dst  = data + slot * MAX_CHANNELS;

        for (int i = 0; i < n_ch; ++i)
            dst[i] = frame.amplitudes[channels[i]];

        // Memory barrier: data must be written before write_pos advances.
        // Python uses write_pos as the "safe to read up to here" marker.
        __sync_synchronize();
        hdr->write_pos    += 1;
        hdr->frame_number  = frame.frameInfo.frame_number;

        ++frames_written;
        // Print heartbeat every second
        if (frames_written % 20000 == 0)
            std::cout << "[streamer] " << frames_written / 20000
                      << " s  frame=" << frame.frameInfo.frame_number
                      << "\r" << std::flush;
    }

    // ── Cleanup ───────────────────────────────────────────────────────────
    maxlab::verifyStatus(maxlab::DataStreamerRaw_close());
    munmap(shm_ptr, SHM_TOTAL);
    shm_unlink(SHM_NAME);
    std::cout << "\n[streamer] Stopped cleanly.\n";
    return 0;
}
