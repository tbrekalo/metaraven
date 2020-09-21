// author tbrekalo 2020

#ifndef RAVEN_COMMON_HPP_
#define RAVEN_COMMON_HPP_

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"

namespace raven {
namespace constants {

char constexpr kFillerSeqsPath[] = "extracted.fasta";

std::size_t constexpr kMinSequenceLen = 1000;

float constexpr kKMerDiscardFreqHard = 0.001;

float constexpr kMerDiscardFreqSoft = 0.00001;

// TODO: Relation between kNonChericLowLimit
std::size_t constexpr kTrimLim = 800;

std::size_t constexpr kSeqsBatchLim = 1ULL << 32;

std::size_t constexpr kOvlpBatchLim = 1ULL << 30;

std::size_t constexpr kFillerLenLim = 20000;

// number of overlaps used per lhs in greedy construction/assembly
std::size_t constexpr kMaxGreedyOvlp = 8; // TODO: consider removing

}  // namespace constants

namespace util {

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path);

// takes ownership of both containers
std::vector<std::unique_ptr<biosoup::Sequence>> MergeSequences(
    std::vector<std::unique_ptr<biosoup::Sequence>>& seqs_a,
    std::vector<std::unique_ptr<biosoup::Sequence>>& seqs_b);

std::vector<std::unique_ptr<biosoup::Sequence>>& NormalizeSeqIds(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences);

std::vector<std::unique_ptr<biosoup::Sequence>> LoadSequences(
    std::string const& path);

std::vector<std::unique_ptr<biosoup::Sequence>> LoadFillerSeqs();

// modifies the original collection
std::vector<std::unique_ptr<biosoup::Sequence>>& TrimSequences(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences);

}  // namespace util

}  // namespace raven

#endif /* RAVEN_COMMON_HPP_ */