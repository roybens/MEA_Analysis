/**
 * mxw_streamer.cpp
 * ================
 * Reads raw frames from MaxOne via DataStreamerRaw_receiveNextFrame,
 * extracts a small set of channels, and writes them into a POSIX
 * shared memory ring buffer that the Python oscilloscope reads.
 *
 * Compile:
 *   g++ -std=c++17 -O3 -o mxw_streamer mxw_streamer.cpp \
 *       -I/path/to/maxlab_lib/include \
 *       -L/path/to/maxlab_lib/lib -lmaxlab -lrt -lpthread \
 *       -Wl,-rpath,/path/to/maxlab_lib/lib
 *
 * Run (AFTER oscilloscope.py has started and created the shared memory):
 *   ./mxw_streamer 42 55 67 88
 *   (pass the channel INDICES — not electrode indices — to stream)
 *
 * The Python oscilloscope.py will tell you the channel indices once
 * it has mapped electrode → channel from the .h5 or config file.
 */

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <thread>
#include <chrono>
#include <atomic>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <semaphore.h>

#include "maxlab/maxlab.h"

// ─────────────────────────────────────────────────────────────────────────────
// Shared memory layout (must match oscilloscope.py exactly)
// ─────────────────────────────────────────────────────────────────────────────

static constexpr const char* SHM_NAME      = "/mxw_scope";
static constexpr int         MAX_CHANNELS  = 8;
static constexpr int         RING_FRAMES   = 40000;  // 2 s @ 20 kHz

// Header sits at offset 0 in shared memory
struct ShmHeader {
    uint32_t magic;           // 0xMEA1CAFE  — Python checks this
    uint32_t n_channels;      // how many channels are active
    uint16_t channel_ids[MAX_CHANNELS];   // which channel indices
    uint64_t write_pos;       // next frame slot to write (wraps mod RING_FRAMES)
    uint64_t frame_number;    // latest hardware frame number written
    uint32_t sample_rate;     // Hz
    float    lsb_uv;          // µV per LSB (default 2.34 for gain=512)
    // 1-bit stim flag per frame: stim happened near this frame
    // Stored as uint8 array, 1 = stim pulse fired at or near this frame
};

static constexpr size_t HEADER_SIZE  = sizeof(ShmHeader);
static constexpr size_t DATA_OFFSET  = 4096;   // data starts here (page-aligned)
// Data layout: RING_FRAMES × MAX_CHANNELS × float32
static constexpr size_t DATA_SIZE    = RING_FRAMES * MAX_CHANNELS * sizeof(float);
static constexpr size_t STIM_OFFSET  = DATA_OFFSET + DATA_SIZE;
static constexpr size_t STIM_SIZE    = RING_FRAMES * sizeof(uint8_t);
static constexpr size_t SHM_TOTAL    = STIM_OFFSET + STIM_SIZE;

static constexpr uint32_t SHM_MAGIC  = 0xC0FFEE42;

static constexpr uint8_t  TARGET_WELL = 0;

// ─────────────────────────────────────────────────────────────────────────────

static std::atomic<bool> g_running{true};

void signal_handler(int) { g_running = false; }

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <ch_idx1> [ch_idx2] [ch_idx3] [ch_idx4]\n"
                  << "  Max " << MAX_CHANNELS << " channel indices.\n"
                  << "  Get indices from oscilloscope.py startup output.\n";
        return 1;
    }

    // Parse channel indices
    std::vector<int> channels;
    for (int i = 1; i < argc && i <= MAX_CHANNELS; ++i)
        channels.push_back(std::atoi(argv[i]));
    int n_ch = (int)channels.size();

    std::cout << "[streamer] Channels: ";
    for (int c : channels) std::cout << c << " ";
    std::cout << "\n";

    // ── Create / open shared memory ───────────────────────────────────────
    shm_unlink(SHM_NAME);   // clean up any stale segment
    int shm_fd = shm_open(SHM_NAME, O_CREAT | O_RDWR, 0666);
    if (shm_fd < 0) { perror("shm_open"); return 1; }
    if (ftruncate(shm_fd, SHM_TOTAL) < 0) { perror("ftruncate"); return 1; }

    void* shm_ptr = mmap(nullptr, SHM_TOTAL, PROT_READ | PROT_WRITE,
                          MAP_SHARED, shm_fd, 0);
    if (shm_ptr == MAP_FAILED) { perror("mmap"); return 1; }
    close(shm_fd);
    memset(shm_ptr, 0, SHM_TOTAL);

    // Fill header
    ShmHeader* hdr = reinterpret_cast<ShmHeader*>(shm_ptr);
    hdr->n_channels  = (uint32_t)n_ch;
    hdr->sample_rate = 20000;
    hdr->lsb_uv      = 2.34f;   // approx µV/LSB at gain 512
    hdr->write_pos   = 0;
    hdr->frame_number = 0;
    for (int i = 0; i < n_ch; ++i)
        hdr->channel_ids[i] = (uint16_t)channels[i];

    float*   data_buf  = reinterpret_cast<float*>  ((char*)shm_ptr + DATA_OFFSET);
    uint8_t* stim_buf  = reinterpret_cast<uint8_t*>((char*)shm_ptr + STIM_OFFSET);

    // Signal magic LAST so Python knows layout is ready
    hdr->magic = SHM_MAGIC;

    std::cout << "[streamer] Shared memory ready: " << SHM_NAME
              << "  (" << (SHM_TOTAL / 1024) << " KB)\n";

    // ── Open MaxOne raw stream ────────────────────────────────────────────
    maxlab::checkVersions();
    maxlab::verifyStatus(maxlab::DataStreamerRaw_open());
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    std::cout << "[streamer] Raw stream open. Streaming...\n";

    signal(SIGINT,  signal_handler);
    signal(SIGTERM, signal_handler);

    maxlab::RawFrameData frame;
    uint64_t frames_written = 0;

    while (g_running) {
        maxlab::Status st = maxlab::DataStreamerRaw_receiveNextFrame(&frame);
        if (st == maxlab::Status::MAXLAB_NO_FRAME) continue;
        if (frame.frameInfo.well_id != TARGET_WELL) continue;
        if (frame.frameInfo.corrupted) continue;

        uint64_t slot = hdr->write_pos % RING_FRAMES;
        float*   dst  = data_buf + slot * MAX_CHANNELS;

        for (int i = 0; i < n_ch; ++i)
            dst[i] = frame.amplitudes[channels[i]];

        stim_buf[slot] = 0;   // Python sets stim flags from its own pulse timing

        // Atomic-ish update: write data BEFORE advancing write_pos
        // (Python reads write_pos to know where the latest data is)
        __sync_synchronize();
        hdr->write_pos    = hdr->write_pos + 1;
        hdr->frame_number = frame.frameInfo.frame_number;

        ++frames_written;
        if (frames_written % 20000 == 0)
            std::cout << "[streamer] " << frames_written / 20000 << " s streamed\r" << std::flush;
    }

    maxlab::verifyStatus(maxlab::DataStreamerRaw_close());
    munmap(shm_ptr, SHM_TOTAL);
    shm_unlink(SHM_NAME);
    std::cout << "\n[streamer] Stopped.\n";
    return 0;
}
