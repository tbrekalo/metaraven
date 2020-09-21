#include <iostream>

#include "common.hpp"

namespace raven {

namespace util {

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [](const std::string& s, const std::string& suff) {
    return s.size() < suff.size()
               ? false
               : s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastaParser>(path);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<
          bioparser::FastqParser>(path);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[raven::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)"
            << std::endl;
  return nullptr;
}

std::vector<std::unique_ptr<biosoup::Sequence>> MergeSequences(
    std::vector<std::unique_ptr<biosoup::Sequence>>& seqs_a,
    std::vector<std::unique_ptr<biosoup::Sequence>>& seqs_b) {
  std::vector<std::unique_ptr<biosoup::Sequence>> dst;

  std::move(seqs_a.begin(), seqs_a.end(), std::back_inserter(dst));
  std::move(seqs_b.begin(), seqs_b.end(), std::back_inserter(dst));

  seqs_a.clear();
  seqs_b.clear();

  return dst;
}

std::vector<std::unique_ptr<biosoup::Sequence>>& NormalizeSeqIds(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {
  biosoup::Sequence::num_objects = sequences.size();
  for (std::size_t i = 0; i < sequences.size(); ++i) {
    sequences[i]->id = i;
  }

  return sequences;
}

std::vector<std::unique_ptr<biosoup::Sequence>> LoadSequences(
    std::string const& path) {
  auto sparser = CreateParser(path);
  auto sequences = sparser->Parse(-1);

  if (sequences.empty()) {
    throw std::runtime_error("[raven::] error: empty sequences set");
  }

  return sequences;
}

std::vector<std::unique_ptr<biosoup::Sequence>> LoadFillerSeqs() {
  return LoadSequences(constants::kFillerSeqsPath);
}

std::vector<std::unique_ptr<biosoup::Sequence>>& TrimSequences(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {
  auto const trim_str = [](std::string& s) -> std::string {
    auto const pos = constants::kTrimLim;
    auto const len = s.size() - 2 * constants::kTrimLim;

    return s.substr(pos, len);
  };

  for (auto& it : sequences) {
    if (it->data.size() > constants::kMinSequenceLen) {
      it->data = trim_str(it->data);
    }
  }

  return sequences;
}

}  // namespace util

}  // namespace raven