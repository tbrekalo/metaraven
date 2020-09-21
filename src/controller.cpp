
#include <iostream>
#include <iterator>
#include <getopt.h>
#include <fstream>

#include "controller.hpp"
#include "common.hpp"

namespace raven {

namespace detail {

static struct option options[] = {
    {"weaken", no_argument, nullptr, 'w'},
    {"polishing-rounds", required_argument, nullptr, 'p'},
    {"match", required_argument, nullptr, 'm'},
    {"mismatch", required_argument, nullptr, 'n'},
    {"gap", required_argument, nullptr, 'g'},
#ifdef CUDA_ENABLED
    {"cuda-poa-batches", optional_argument, nullptr, 'c'},
    {"cuda-banded-alignment", no_argument, nullptr, 'b'},
    {"cuda-alignment-batches", required_argument, nullptr, 'a'},
#endif
    {"graphical-fragment-assembly", required_argument, nullptr, 'f'},
    {"second-run", no_argument, nullptr, 's'},
    {"resume", no_argument, nullptr, 'r'},
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}};

void RavenRun(Config const& conf, Data& data) {
  data.graph.Construct(data.sequences);
  data.graph.Assemble();
  data.graph.Polish(data.sequences, conf.m, conf.n, conf.g,
                    conf.cuda_poa_batches, conf.cuda_banded_alignment,
                    conf.cuda_alignment_batches, conf.num_polishing_rounds);
}

void RavenPrintResults(Config const& conf, Data& data) {
  data.graph.PrintGFA(conf.gfa_path);

  for (auto const& it : data.graph.GetUnitigs(conf.num_polishing_rounds > 0)) {
    std::cout << ">" << it->name << std::endl;
    std::cout << it->data << std::endl;
  }

  data.timer.Stop();
  std::cerr << "[raven::] " << std::fixed << data.timer.elapsed_time() << "s"
            << std::endl;
}

void SingleRun(Config const& conf, Data& data) {
  RavenRun(conf, data);
  RavenPrintResults(conf, data);
}

void DoubleRun(Config const& conf, Data& data) {
  RavenRun(conf, data);

  std::cerr << "[raven::] finished first run in: " << std::fixed
            << data.timer.Stop() << "s\n"
            << std::endl;

  std::vector<std::unique_ptr<biosoup::Sequence>> unitigs;
  for (auto const& it : data.graph.GetUnitigs(conf.num_polishing_rounds > 0)) {
    unitigs.emplace_back(new biosoup::Sequence(">" + it->name, it->data));
  }

  data.graph.Clear();
  data.timer.Start();

  auto const expected = data.graph.GreedyConstruct(util::NormalizeSeqIds(unitigs));
  data.graph.GreedyAssemble(expected);


  RavenPrintResults(conf, data);
}

}  // namespace detail

Config ParseConfig(int argc, char** argv) {
  Config conf;
  std::string optstr = "p:m:n:g:t:h";

#ifdef CUDA_ENABLED
  optstr += "c:ba:";
#endif
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), detail::options,
                            nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 'w':
        conf.weaken = true;
        break;
      case 'p':
        conf.num_polishing_rounds = atoi(optarg);
        break;
      case 'm':
        conf.m = atoi(optarg);
        break;
      case 'n':
        conf.n = atoi(optarg);
        break;
      case 'g':
        conf.g = atoi(optarg);
        break;
#ifdef CUDA_ENABLED
      case 'c':
        conf.cuda_poa_batches = 1;
        // next text entry is not an option, assuming it's the arg for option c
        if (optarg == NULL && argv[optind] != NULL && argv[optind][0] != '-') {
          conf.cuda_poa_batches = atoi(argv[optind++]);
        }
        // optional argument provided in the ususal way
        if (optarg != NULL) {
          conf.cuda_poa_batches = atoi(optarg);
        }
        break;
      case 'b':
        conf.cuda_banded_alignment = true;
        break;
      case 'a':
        conf.cuda_alignment_batches = atoi(optarg);
        break;
#endif
      case 'f':
        conf.gfa_path = optarg;
        break;
      case 'r':
        conf.resume = true;
        break;
      case 's':
        conf.second_run = true;
        break;
      case 't':
        conf.num_threads = atoi(optarg);
        break;
      case 'v':
        std::cout << raven_version << std::endl;
        conf.run = false;
        break;
      case 'h':
        Help();
        conf.run = false;
      default:
        break;
    }
  }

  if (conf.run && optind >= argc) {
    Help();
    throw std::runtime_error("[raven::] error: missing input file!");
  } else if (conf.run) {
    conf.sequence_path = argv[optind];
  }

  return conf;
}

void Help() {
  std::cout
      << "usage: metaraven [options ...] <sequences>\n"
         "\n"
         "  # default output is stdout\n"
         "  <sequences>\n"
         "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
         "\n"
         "  options:\n"
         "    --weaken\n"
         "      use larger (k, w) when assembling highly accurate sequences\n"
         "    -p, --polishing-rounds <int>\n"
         "      default: 2\n"
         "      number of times racon is invoked\n"
         "    -m, --match <int>\n"
         "      default: 3\n"
         "      score for matching bases\n"
         "    -n, --mismatch <int>\n"
         "      default: -5\n"
         "      score for mismatching bases\n"
         "    -g, --gap <int>\n"
         "      default: -4\n"
         "      gap penalty (must be negative)\n"
#ifdef CUDA_ENABLED
         "    -c, --cuda-poa-batches <int>\n"
         "       default: 0\n"
         "       number of batches for CUDA accelerated polishing\n"
         "    -b, --cuda-banded-alignment\n"
         "       use banding approximation for polishing on GPU\n"
         "       (only applicable when -c is used)\n"
         "    -a, --cuda-alignment-batches <int>\n"
         "       default: 0\n"
         "       number of batches for CUDA accelerated alignment\n"
#endif
         "    --second-run\n"
         "      reuses non-chimeric in combination with unitigs"
         "    --graphical-fragment-assembly <string>\n"
         "      prints the assemblg graph in GFA format\n"
         "    --resume\n"
         "      resume previous run from last checkpoint\n"
         "    -t, --threads <int>\n"
         "      default: 1\n"
         "      number of threads\n"
         "    --version\n"
         "      prints the version number\n"
         "    -h, --help\n"
         "       prints the usage\n";
}

Data::Data(Config const& conf)
    : sequences{},
      thread_pool{std::make_shared<thread_pool::ThreadPool>(conf.num_threads)},
      graph{conf.weaken, thread_pool} {
  timer.Start();
}

Data::Data(bool const weaken, std::uint32_t const num_threads)
    : sequences{},
      thread_pool{std::make_shared<thread_pool::ThreadPool>(num_threads)},
      graph{weaken, thread_pool} {
  timer.Start();
}

Data Setup(Config const& conf) {
  Data data(conf);
  auto sparser = util::CreateParser(conf.sequence_path);

  if (conf.resume) {
    data.graph.Load();

    std::cerr << "[raven::] loaded previous run " << std::fixed
              << data.timer.Stop() << "s" << std::endl;
  }

  if (conf.second_run) {
    std::ofstream os{constants::kFillerSeqsPath, std::ios_base::trunc};
  }

  if (data.graph.stage() < -3 ||
      conf.num_polishing_rounds > std::max(0, data.graph.stage())) {
    data.sequences = util::LoadSequences(conf.sequence_path);

    std::cerr << "[raven::] loaded " << data.sequences.size() << " sequences "
              << std::fixed << data.timer.Stop() << "s" << std::endl;

    data.timer.Start();
  }

  return data;
}  // namespace raven

void RunRaven(Config const& conf, Data& data) {
  if (!conf.second_run) {
    detail::SingleRun(conf, data);
  } else {
    detail::DoubleRun(conf, data);
  }
}

}  // namespace raven