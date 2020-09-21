// Author tbrekalo 2020

#ifndef RAVEN_CONTORLLER_HPP_
#define RAVEN_CONTROLLER_HPP_

#include <stdexcept>
#include <cstdint>
#include <string>
#include <memory>
#include <vector>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"
#include "biosoup/timer.hpp"

#include "graph.hpp"

namespace raven {

static const char* raven_version = RAVEN_VERSION;

struct Config {
  bool run = true;
  bool second_run = false;

  bool weaken = false;

  std::int32_t num_polishing_rounds = 2;
  std::int8_t m = 3;
  std::int8_t n = -5;
  std::int8_t g = -4;

  std::string gfa_path = "";
  bool resume = false;

  std::uint32_t num_threads = 1;

  std::uint32_t cuda_poa_batches = 0;
  std::uint32_t cuda_alignment_batches = 0;
  bool cuda_banded_alignment = false;

  std::string sequence_path;
};

struct Data {
public:
  Data(Config const& conf);
  Data(bool weaken, std::uint32_t num_threads);

  std::vector<std::unique_ptr<biosoup::Sequence>> sequences;
  std::shared_ptr<thread_pool::ThreadPool> thread_pool;
  raven::Graph graph;

  biosoup::Timer timer{};
};

Config ParseConfig(int argc, char** argv);

void Help();

Data Setup(Config const&);

void RunRaven(Config const&, Data&);

}  // namespace raven

#endif  // RAVEN_CONTOLLER_HPP_