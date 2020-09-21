// Copyright (c) 2020 Robert Vaser

#include <iostream>
#include <cstdlib>

#include <getopt.h>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/timer.hpp"

#include "controller.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

int main(int argc, char** argv) {
  try {
    auto conf = raven::ParseConfig(argc, argv);
    if (conf.run) {
      auto data = raven::Setup(conf);
      raven::RunRaven(conf, data);
    }

  } catch (std::exception const& e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
