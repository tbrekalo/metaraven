// Copyright (c) 2020 Robert Vaser

#include "graph.hpp"

#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <exception>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <numeric>
#include <random>
#include <deque>

#include <assert.h>  // TODO: Remove after dev

#include "cereal/archives/binary.hpp"
#include "cereal/archives/json.hpp"
#include "racon/polisher.hpp"
#include "biosoup/timer.hpp"

#include "common.hpp"

namespace raven {

namespace detail {

biosoup::Overlap OverlapReverse(biosoup::Overlap const& o) {
  return biosoup::Overlap(o.rhs_id, o.rhs_begin, o.rhs_end, o.lhs_id,
                          o.lhs_begin, o.lhs_end, o.score, o.strand);
}

std::uint32_t OverlapLength(biosoup::Overlap const& o) {
  return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
}

bool OverlapUpdate(std::vector<std::unique_ptr<Pile>>& piles,
                   biosoup::Overlap& o) {
  if (piles[o.lhs_id]->is_invalid() || piles[o.rhs_id]->is_invalid()) {
    return false;
  }
  if (o.lhs_begin >= piles[o.lhs_id]->end() ||
      o.lhs_end <= piles[o.lhs_id]->begin() ||
      o.rhs_begin >= piles[o.rhs_id]->end() ||
      o.rhs_end <= piles[o.rhs_id]->begin()) {
    return false;
  }

  std::uint32_t lhs_begin =
      o.lhs_begin + (o.strand ? (o.rhs_begin < piles[o.rhs_id]->begin()
                                     ? piles[o.rhs_id]->begin() - o.rhs_begin
                                     : 0)
                              : (o.rhs_end > piles[o.rhs_id]->end()
                                     ? o.rhs_end - piles[o.rhs_id]->end()
                                     : 0));
  std::uint32_t lhs_end =
      o.lhs_end - (o.strand ? (o.rhs_end > piles[o.rhs_id]->end()
                                   ? o.rhs_end - piles[o.rhs_id]->end()
                                   : 0)
                            : (o.rhs_begin < piles[o.rhs_id]->begin()
                                   ? piles[o.rhs_id]->begin() - o.rhs_begin
                                   : 0));

  std::uint32_t rhs_begin =
      o.rhs_begin + (o.strand ? (o.lhs_begin < piles[o.lhs_id]->begin()
                                     ? piles[o.lhs_id]->begin() - o.lhs_begin
                                     : 0)
                              : (o.lhs_end > piles[o.lhs_id]->end()
                                     ? o.lhs_end - piles[o.lhs_id]->end()
                                     : 0));
  std::uint32_t rhs_end =
      o.rhs_end - (o.strand ? (o.lhs_end > piles[o.lhs_id]->end()
                                   ? o.lhs_end - piles[o.lhs_id]->end()
                                   : 0)
                            : (o.lhs_begin < piles[o.lhs_id]->begin()
                                   ? piles[o.lhs_id]->begin() - o.lhs_begin
                                   : 0));

  if (lhs_begin >= piles[o.lhs_id]->end() ||
      lhs_end <= piles[o.lhs_id]->begin() ||
      rhs_begin >= piles[o.rhs_id]->end() ||
      rhs_end <= piles[o.rhs_id]->begin()) {
    return false;
  }

  lhs_begin = std::max(lhs_begin, piles[o.lhs_id]->begin());
  lhs_end = std::min(lhs_end, piles[o.lhs_id]->end());
  rhs_begin = std::max(rhs_begin, piles[o.rhs_id]->begin());
  rhs_end = std::min(rhs_end, piles[o.rhs_id]->end());

  if (lhs_begin >= lhs_end || lhs_end - lhs_begin < 84 ||
      rhs_begin >= rhs_end || rhs_end - rhs_begin < 84) {
    return false;
  }

  o.lhs_begin = lhs_begin;
  o.lhs_end = lhs_end;
  o.rhs_begin = rhs_begin;
  o.rhs_end = rhs_end;

  return true;
}

std::uint32_t OverlapType(std::vector<std::unique_ptr<Pile>>& piles,
                          biosoup::Overlap const& o) {
  std::uint32_t lhs_length = piles[o.lhs_id]->end() - piles[o.lhs_id]->begin();
  std::uint32_t lhs_begin = o.lhs_begin - piles[o.lhs_id]->begin();
  std::uint32_t lhs_end = o.lhs_end - piles[o.lhs_id]->begin();

  std::uint32_t rhs_length = piles[o.rhs_id]->end() - piles[o.rhs_id]->begin();
  std::uint32_t rhs_begin =
      o.strand ? o.rhs_begin - piles[o.rhs_id]->begin()
               : rhs_length - (o.rhs_end - piles[o.rhs_id]->begin());
  std::uint32_t rhs_end =
      o.strand ? o.rhs_end - piles[o.rhs_id]->begin()
               : rhs_length - (o.rhs_begin - piles[o.rhs_id]->begin());

  std::uint32_t overhang = std::min(lhs_begin, rhs_begin) +
                           std::min(lhs_length - lhs_end, rhs_length - rhs_end);

  if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
      rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
    return 0;  // internal
  }
  if (lhs_begin <= rhs_begin && lhs_length - lhs_end <= rhs_length - rhs_end) {
    return 1;  // lhs contained
  }
  if (rhs_begin <= lhs_begin && rhs_length - rhs_end <= lhs_length - lhs_end) {
    return 2;  // rhs contained
  }
  if (lhs_begin > rhs_begin) {
    return 3;  // lhs -> rhs
  }
  return 4;  // rhs -> lhs
}

bool OverlapFinalize(std::vector<std::unique_ptr<Pile>>& piles,
                     biosoup::Overlap& o) {
  o.score = OverlapType(piles, o);
  if (o.score < 3) {
    return false;
  }

  o.lhs_begin -= piles[o.lhs_id]->begin();
  o.lhs_end -= piles[o.lhs_id]->begin();

  o.rhs_begin -= piles[o.rhs_id]->begin();
  o.rhs_end -= piles[o.rhs_id]->begin();
  if (!o.strand) {
    auto rhs_begin = o.rhs_begin;
    o.rhs_begin = piles[o.rhs_id]->length() - o.rhs_end;
    o.rhs_end = piles[o.rhs_id]->length() - rhs_begin;
  }
  return true;
}

// prerequsite, sequence id corresponds to pile id
void StoreValidRegions(
    std::vector<std::unique_ptr<Pile>>& piles,
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {
  std::size_t cnt = 0;
  std::ofstream os(constants::kFillerSeqsPath);

  for (auto const& seq : sequences) {
    auto const& pile = piles[seq->id];

    auto const pos = pile->begin();
    auto const len = pile->length();

    if (len >= constants::kMinSequenceLen) {
      auto const valid_subsequence = seq->data.substr(pos, len);
      ++cnt;

      os << ">nc" + std::to_string(seq->id) << '\n'
         << valid_subsequence << '\n';
    }
  }

  std::cerr << "[raven::Graph::detail::StoreValidRegions] saved " << cnt
            << " sequence regions" << std::endl;
}

enum class OverlapSide : std::uint8_t { kLeft, kRight };

// assumes that sequence id corresponds to ovlp.lhs_id
OverlapSide DetermineOverlapSide(
    std::unique_ptr<biosoup::Sequence> const& sequence,
    biosoup::Overlap const& ovlp) {
  if (ovlp.lhs_begin < sequence->data.size() / 2) {
    return OverlapSide::kLeft;
  } else {
    return OverlapSide::kRight;
  }
}

class OverlapSideVecs {
 public:
  void AddOverlap(std::unique_ptr<biosoup::Sequence> const& sequence,
                  biosoup::Overlap const& ovlp) {
    auto const side = DetermineOverlapSide(sequence, ovlp);
    if (side == OverlapSide::kLeft) {
      AddTo(left_, ovlp);
    } else {
      AddTo(right_, ovlp);
    }
  }

  // moves container content into a new vector
  std::vector<biosoup::Overlap> MergedSides() {
    std::vector<biosoup::Overlap> dst;

    std::move(left_.begin(), left_.end(), std::back_inserter(dst));
    std::move(right_.begin(), right_.end(), std::back_inserter(dst));

    return dst;
  }

 private:
  void AddTo(std::vector<biosoup::Overlap>& vec, biosoup::Overlap const& ovlp) {
    // assumes small vector size
    auto const iter_pos = std::find_if(
        vec.begin(), vec.end(), [&ovlp](biosoup::Overlap const& val) -> bool {
          return OverlapLength(val) < OverlapLength(ovlp);
        });

    vec.emplace(iter_pos, ovlp);
    if (vec.size() > constants::kMaxGreedyOvlp) {
      vec.resize(constants::kMaxGreedyOvlp);
    }
  }

  std::vector<biosoup::Overlap> left_;
  std::vector<biosoup::Overlap> right_;
};

enum class OverlapCategory { kIrrelevant, kLeft, kRight };

enum class ExpandDir { kLeft, kRight };

}  // namespace detail

Graph::Node::Node(const biosoup::Sequence& sequence)
    : id(num_objects++),
      name(sequence.name),
      data(sequence.data),
      count(1),
      is_circular(),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {}

Graph::Node::Node(Node* begin, Node* end)
    : id(num_objects++),
      name(),
      data(),
      count(),
      is_circular(begin == end),
      is_polished(),
      transitive(),
      inedges(),
      outedges(),
      pair() {
  auto it = begin;
  while (true) {
    data += it->outedges.front()->Label();
    count += it->count;
    if ((it = it->outedges.front()->head) == end) {
      break;
    }
  }
  if (begin != end) {
    data += end->data;
    count += end->count;
  }

  name = (is_unitig() ? "Utg" : "Ctg") + std::to_string(id);
}

Graph::Edge::Edge(Node* tail, Node* head, std::uint32_t length)
    : id(num_objects++),
      length(length),
      weight(0),
      tail(tail),
      head(head),
      pair() {
  tail->outedges.emplace_back(this);
  head->inedges.emplace_back(this);
}

std::atomic<std::uint32_t> Graph::Node::num_objects{0};
std::atomic<std::uint32_t> Graph::Edge::num_objects{0};

Graph::Graph(bool weaken, std::shared_ptr<thread_pool::ThreadPool> thread_pool)
    : thread_pool_(thread_pool ? thread_pool
                               : std::make_shared<thread_pool::ThreadPool>(1)),
      minimizer_engine_(weaken ? 29 : 15, weaken ? 9 : 5, thread_pool_),
      stage_(-5),
      piles_(),
      nodes_(),
      edges_() {}

std::vector<std::unique_ptr<biosoup::Sequence>> Graph::Preprocess(
    std::vector<std::unique_ptr<biosoup::Sequence>>&& sequences) {
  // container for new sequences
  std::vector<std::unique_ptr<biosoup::Sequence>> dst;

  biosoup::Timer timer{};
  std::vector<std::vector<biosoup::Overlap>> overlaps{};

  timer.Start();
  std::size_t sequence_batch_bytes = 0;
  for (std::size_t j = 0, i = 0; i < sequences.size(); ++i) {
    sequence_batch_bytes += sequences[i]->data.size();
    if (i + 1 != sequences.size() &&
        sequence_batch_bytes < constants::kSeqsBatchLim) {
      continue;
    }

    minimizer_engine_.Minimize(sequences.begin() + j,
                               sequences.begin() + i + 1);
    minimizer_engine_.Filter(constants::kKMerDiscardFreqHard);

    std::cerr << "[raven::Graph::Preprocess] minimized " << j << " - " << i + 1
              << " / " << sequences.size() << " " << std::fixed << timer.Stop()
              << "s" << std::endl;

    std::vector<std::future<std::vector<biosoup::Overlap>>> overlap_futures;

    std::size_t overlap_batch_bytes = 0;
    for (std::size_t k = j; k < i + 1; ++k) {
      overlap_futures.emplace_back(thread_pool_->Submit(
          [&](std::size_t pos) -> std::vector<biosoup::Overlap> {
            return minimizer_engine_.Map(sequences[pos], true, true, false);
          },
          k));

      overlap_batch_bytes += sequences[i]->data.size();
      if (k != i && overlap_batch_bytes < constants::kOvlpBatchLim) {
        continue;
      }

      for (auto& future : overlap_futures) {
        for (auto& ovlp : future.get()) {
          overlaps[ovlp.lhs_id].emplace_back(ovlp);
          overlaps[ovlp.rhs_id].emplace_back(detail::OverlapReverse(ovlp));
        }
      }

      overlap_futures.clear();

      j = i + 1;
    }
  }
  Clear();  // reset information before main algorithm
  return dst;
}

void Graph::Construct(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {  // NOLINT

  if (sequences.empty() || stage_ > -4) {
    return;
  }

  std::vector<std::vector<biosoup::Overlap>> overlaps(sequences.size());

  // biosoup::Overlap helper functions
  auto overlap_reverse = [](const biosoup::Overlap& o) -> biosoup::Overlap {
    return biosoup::Overlap(o.rhs_id, o.rhs_begin, o.rhs_end, o.lhs_id,
                            o.lhs_begin, o.lhs_end, o.score, o.strand);
  };
  auto overlap_length = [](const biosoup::Overlap& o) -> std::uint32_t {
    return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
  };

  auto overlap_update = [&](biosoup::Overlap& o) -> bool {
    return detail::OverlapUpdate(piles_, o);
  };

  auto overlap_type = [&](const biosoup::Overlap& o) -> std::uint32_t {
    return detail::OverlapType(piles_, o);
  };

  auto overlap_finalize = [&](biosoup::Overlap& o) -> bool {
    return detail::OverlapFinalize(piles_, o);
  };

  auto connected_components =
      [&]() -> std::vector<std::vector<std::uint32_t>> {  // NOLINT
    std::vector<std::vector<std::uint32_t>> connections(sequences.size());
    for (const auto& it : overlaps) {
      for (const auto& jt : it) {
        if (overlap_type(jt) > 2) {
          connections[jt.lhs_id].emplace_back(jt.rhs_id);
          connections[jt.rhs_id].emplace_back(jt.lhs_id);
        }
      }
    }

    std::vector<std::vector<std::uint32_t>> dst;
    std::vector<char> is_visited(sequences.size(), false);
    for (std::uint32_t i = 0; i < connections.size(); ++i) {
      if (piles_[i]->is_invalid() || is_visited[i]) {
        continue;
      }

      dst.resize(dst.size() + 1);
      std::deque<std::uint32_t> que = {i};
      while (!que.empty()) {
        std::uint32_t j = que.front();
        que.pop_front();

        if (is_visited[j]) {
          continue;
        }
        is_visited[j] = true;
        dst.back().emplace_back(j);

        for (const auto& it : connections[j]) {
          que.emplace_back(it);
        }
      }
    }

    return dst;
  };
  // biosoup::Overlap helper functions

  if (stage_ == -5) {  // checkpoint test
    Store();
  }

  biosoup::Timer timer{};

  if (stage_ == -5) {  // find overlaps and create piles
    for (const auto& it : sequences) {
      piles_.emplace_back(new Pile(it->id, it->data.size()));
    }
    std::size_t bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < sequences.size(); ++i) {
      bytes += sequences[i]->data.size();
      if (i != sequences.size() - 1 && bytes < constants::kSeqsBatchLim) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(sequences.begin() + j,
                                 sequences.begin() + i + 1, true);
      minimizer_engine_.Filter(constants::kKMerDiscardFreqHard);

      std::cerr << "[raven::Graph::Construct] minimized " << j << " - " << i + 1
                << " / " << sequences.size() << " " << std::fixed
                << timer.Stop() << "s" << std::endl;

      timer.Start();

      std::vector<std::uint32_t> num_overlaps(overlaps.size());
      for (std::uint32_t k = 0; k < overlaps.size(); ++k) {
        num_overlaps[k] = overlaps[k].size();
      }

      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;

      for (std::uint32_t k = 0; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, true, true);
            },
            k));

        bytes += sequences[k]->data.size();
        if (k != i && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto& it : thread_futures) {
          for (const auto& jt : it.get()) {
            overlaps[jt.lhs_id].emplace_back(jt);
            overlaps[jt.rhs_id].emplace_back(overlap_reverse(jt));
          }
        }
        thread_futures.clear();

        std::vector<std::future<void>> void_futures;
        for (const auto& it : piles_) {
          if (overlaps[it->id()].empty() ||
              overlaps[it->id()].size() == num_overlaps[it->id()]) {
            continue;
          }

          void_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                piles_[i]->AddLayers(overlaps[i].begin() + num_overlaps[i],
                                     overlaps[i].end());

                num_overlaps[i] =
                    std::min(overlaps[i].size(), static_cast<std::size_t>(16));

                if (overlaps[i].size() < 16) {
                  return;
                }

                std::sort(overlaps[i].begin(), overlaps[i].end(),
                          [&](const biosoup::Overlap& lhs,
                              const biosoup::Overlap& rhs) -> bool {
                            return overlap_length(lhs) > overlap_length(rhs);
                          });

                std::vector<biosoup::Overlap> tmp;
                tmp.insert(tmp.end(), overlaps[i].begin(),
                           overlaps[i].begin() + 16);  // NOLINT
                tmp.swap(overlaps[i]);
              },
              it->id()));
        }
        for (const auto& it : void_futures) {
          it.wait();
        }
      }

      std::cerr << "[raven::Graph::Construct] mapped sequences " << std::fixed
                << timer.Stop() << "s" << std::endl;

      j = i + 1;
    }
  }

  if (stage_ == -5) {  // trim and annotate piles
    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      thread_futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> void {
            piles_[i]->FindValidRegion(4);
            if (piles_[i]->is_invalid()) {
              std::vector<biosoup::Overlap>().swap(overlaps[i]);
            } else {
              piles_[i]->FindMedian();
              piles_[i]->FindChimericRegions();
            }
          },
          i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] annotated piles " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -5) {  // resolve contained reads
    timer.Start();

    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
      std::uint32_t k = 0;
      for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
        if (!overlap_update(overlaps[i][j])) {
          continue;
        }
        std::uint32_t type = overlap_type(overlaps[i][j]);
        if (type == 1 && !piles_[overlaps[i][j].rhs_id]->is_maybe_chimeric()) {
          piles_[i]->set_is_contained();
        } else if (type == 2 && !piles_[i]->is_maybe_chimeric()) {
          piles_[overlaps[i][j].rhs_id]->set_is_contained();
        } else {
          overlaps[i][k++] = overlaps[i][j];
        }
      }
      overlaps[i].resize(k);
    }
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_contained()) {
        piles_[i]->set_is_invalid();
        std::vector<biosoup::Overlap>().swap(overlaps[i]);
      }
    }

    std::cerr << "[raven::Graph::Construct] removed contained sequences "
              << std::fixed << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -5) {  // resolve chimeric sequences
    timer.Start();

    while (true) {
      auto components = connected_components();
      for (const auto& it : components) {
        std::vector<std::uint32_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                         medians.end());
        std::uint32_t median = medians[medians.size() / 2];

        std::vector<std::future<void>> thread_futures;
        for (const auto& jt : it) {
          thread_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                piles_[i]->ClearChimericRegions(median);
                if (piles_[i]->is_invalid()) {
                  std::vector<biosoup::Overlap>().swap(overlaps[i]);
                }
              },
              jt));
        }
        for (const auto& it : thread_futures) {
          it.wait();
        }
        thread_futures.clear();
      }

      bool is_changed = false;
      for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::uint32_t k = 0;
        for (std::uint32_t j = 0; j < overlaps[i].size(); ++j) {
          if (overlap_update(overlaps[i][j])) {
            overlaps[i][k++] = overlaps[i][j];
          } else {
            is_changed = true;
          }
        }
        overlaps[i].resize(k);
      }

      if (!is_changed) {
        for (const auto& it : overlaps) {
          for (const auto& jt : it) {
            std::uint32_t type = overlap_type(jt);
            if (type == 1) {
              piles_[jt.lhs_id]->set_is_contained();
              piles_[jt.lhs_id]->set_is_invalid();
            } else if (type == 2) {
              piles_[jt.rhs_id]->set_is_contained();
              piles_[jt.rhs_id]->set_is_invalid();
            }
          }
        }
        overlaps.clear();
        break;
      }
    }

    std::cerr << "[raven::Graph::Construct] removed chimeric sequences "
              << std::fixed << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -5) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Construct] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -4) {  // clear piles for sensitive overlaps
    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_invalid()) {
        continue;
      }
      thread_futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> void { piles_[i]->ClearValidRegion(); }, i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] cleared piles " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -4) {  // find overlaps and update piles with repetitive regions
    std::sort(sequences.begin(), sequences.end(),
              [&](const std::unique_ptr<biosoup::Sequence>& lhs,
                  const std::unique_ptr<biosoup::Sequence>& rhs) -> bool {
                return piles_[lhs->id]->is_invalid() <
                           piles_[rhs->id]->is_invalid() ||  // NOLINT
                       (piles_[lhs->id]->is_invalid() ==
                            piles_[rhs->id]->is_invalid() &&
                        lhs->id < rhs->id);  // NOLINT
              });

    std::uint32_t s = 0;
    for (std::uint32_t i = 0; i < sequences.size(); ++i) {
      if (piles_[sequences[i]->id]->is_invalid()) {
        s = i;
        break;
      }
    }

    // map invalid reads to valid reads
    overlaps.resize(sequences.size() + 1);
    std::size_t bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < s; ++i) {
      bytes += sequences[i]->data.size();
      if (i != s - 1 && bytes < (1ULL << 32)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(sequences.begin() + j,
                                 sequences.begin() + i + 1, true);

      std::cerr << "[raven::Graph::Construct] minimized " << j << " - " << i + 1
                << " / " << s << " " << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      minimizer_engine_.Filter(constants::kMerDiscardFreqSoft);
      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
      for (std::uint32_t k = s; k < sequences.size(); ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, false, true);
            },
            k));

        bytes += sequences[k]->data.size();
        if (k != sequences.size() - 1 && bytes < (1U << 30)) {
          continue;
        }
        bytes = 0;

        for (auto& it : thread_futures) {
          for (const auto& jt : it.get()) {
            overlaps[jt.rhs_id].emplace_back(jt);
          }
        }
        thread_futures.clear();

        std::vector<std::future<void>> void_futures;
        for (std::uint32_t k = j; k < i + 1; ++k) {
          if (overlaps[sequences[k]->id].empty()) {
            continue;
          }
          void_futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                piles_[i]->AddLayers(overlaps[i].begin(), overlaps[i].end());
                std::vector<biosoup::Overlap>().swap(overlaps[i]);
              },
              sequences[k]->id));
        }
        for (const auto& it : void_futures) {
          it.wait();
        }
      }

      std::cerr << "[raven::Graph::Construct] mapped invalid sequences "
                << std::fixed << timer.Stop() << "s" << std::endl;

      j = i + 1;
    }

    // map valid reads to each other
    bytes = 0;
    for (std::uint32_t i = 0, j = 0; i < s; ++i) {
      bytes += sequences[i]->data.size();
      if (i != s - 1 && bytes < (1U << 30)) {
        continue;
      }
      bytes = 0;

      timer.Start();

      minimizer_engine_.Minimize(sequences.begin() + j,
                                 sequences.begin() + i + 1);

      std::cerr << "[raven::Graph::Construct] minimized " << j << " - " << i + 1
                << " / " << s << " " << std::fixed << timer.Stop() << "s"
                << std::endl;

      timer.Start();

      std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
      minimizer_engine_.Filter(constants::kKMerDiscardFreqHard);
      for (std::uint32_t k = 0; k < i + 1; ++k) {
        thread_futures.emplace_back(thread_pool_->Submit(
            [&](std::uint32_t i) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[i], true, true);
            },
            k));
      }
      for (auto& it : thread_futures) {
        for (auto& jt : it.get()) {
          if (!overlap_update(jt)) {
            continue;
          }
          std::uint32_t type = overlap_type(jt);
          if (type == 0) {
            continue;
          } else if (type == 1) {
            piles_[jt.lhs_id]->set_is_contained();
          } else if (type == 2) {
            piles_[jt.rhs_id]->set_is_contained();
          } else {
            if (overlaps.back().size() &&
                overlaps.back().back().lhs_id == jt.lhs_id &&
                overlaps.back().back().rhs_id == jt.rhs_id) {
              if (overlap_length(overlaps.back().back()) < overlap_length(jt)) {
                overlaps.back().back() = jt;
              }
            } else {
              overlaps.back().emplace_back(jt);
            }
          }
        }
      }
      thread_futures.clear();

      std::cerr << "[raven::Graph::Construct] mapped valid sequences "
                << std::fixed << timer.Stop() << "s" << std::endl;

      j = i + 1;
    }

    timer.Start();

    std::vector<std::future<void>> thread_futures;
    for (std::uint32_t i = 0; i < piles_.size(); ++i) {
      if (piles_[i]->is_contained()) {
        piles_[i]->set_is_invalid();
        continue;
      }
      if (piles_[i]->is_invalid()) {
        continue;
      }
      thread_futures.emplace_back(thread_pool_->Submit(
          [&](std::uint32_t i) -> void {
            piles_[i]->ClearInvalidRegion();
            piles_[i]->FindMedian();
          },
          i));
    }
    for (const auto& it : thread_futures) {
      it.wait();
    }
    thread_futures.clear();

    std::cerr << "[raven::Graph::Construct] updated piles " << std::fixed
              << timer.Stop() << "s" << std::endl;

    timer.Start();

    {
      std::uint32_t k = 0;
      for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
        if (overlap_update(overlaps.back()[i])) {
          overlaps.back()[k++] = overlaps.back()[i];
        }
      }
      overlaps.back().resize(k);
    }

    std::cerr << "[raven::Graph::Construct] updated overlaps " << std::fixed
              << timer.Stop() << "s" << std::endl;

    std::sort(sequences.begin(), sequences.end(),
              [&](const std::unique_ptr<biosoup::Sequence>& lhs,
                  const std::unique_ptr<biosoup::Sequence>& rhs) -> bool {
                return lhs->id < rhs->id;
              });
  }

  if (stage_ == -4) {  // resolve repeat induced overlaps
    timer.Start();

    while (true) {
      auto components = connected_components();
      for (const auto& it : components) {
        std::vector<std::uint32_t> medians;
        for (const auto& jt : it) {
          medians.emplace_back(piles_[jt]->median());
        }
        std::nth_element(medians.begin(), medians.begin() + medians.size() / 2,
                         medians.end());
        std::uint32_t median = medians[medians.size() / 2];

        std::vector<std::future<void>> futures;
        for (const auto& jt : it) {
          futures.emplace_back(thread_pool_->Submit(
              [&](std::uint32_t i) -> void {
                piles_[i]->FindRepetitiveRegions(median);
              },
              jt));
        }
        for (const auto& it : futures) {
          it.wait();
        }
      }

      for (const auto& it : overlaps.back()) {
        piles_[it.lhs_id]->UpdateRepetitiveRegions(it);
        piles_[it.rhs_id]->UpdateRepetitiveRegions(it);
      }

      bool is_changed = false;
      std::uint32_t j = 0;
      for (std::uint32_t i = 0; i < overlaps.back().size(); ++i) {
        const auto& it = overlaps.back()[i];
        if (piles_[it.lhs_id]->CheckRepetitiveRegions(it) ||
            piles_[it.rhs_id]->CheckRepetitiveRegions(it)) {
          is_changed = true;
        } else {
          overlaps.back()[j++] = it;
        }
      }
      overlaps.back().resize(j);

      if (!is_changed) {
        break;
      }

      for (const auto& it : components) {
        for (const auto& jt : it) {
          piles_[jt]->ClearRepetitiveRegions();
        }
      }
    }

    std::cerr << "[raven::Graph::Construct] removed false overlaps "
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();
  }

  detail::StoreValidRegions(piles_, sequences);

  assert(Node::num_objects == 0);  // TODO: Remove
  if (stage_ == -4) {              // construct assembly graph
    std::vector<std::int32_t> sequence_to_node(piles_.size(), -1);
    for (const auto& it : piles_) {  // create nodes
      if (it->is_invalid()) {
        continue;
      }

      auto sequence = biosoup::Sequence{
          sequences[it->id()]->name,
          sequences[it->id()]->data.substr(it->begin(),
                                           it->end() - it->begin())};  // NOLINT

      sequence_to_node[it->id()] = Node::num_objects;

      auto node = std::make_shared<Node>(sequence);
      sequence.ReverseAndComplement();
      nodes_.emplace_back(node);
      nodes_.emplace_back(std::make_shared<Node>(sequence));
      node->pair = nodes_.back().get();
      node->pair->pair = node.get();
    }

    std::cerr << "[raven::Graph::Construct] stored " << nodes_.size()
              << " nodes "  // NOLINT
              << std::fixed << timer.Stop() << "s" << std::endl;

    timer.Start();

    for (auto& it : overlaps.back()) {  // create edges
      if (!overlap_finalize(it)) {
        continue;
      }

      auto tail = nodes_[sequence_to_node[it.lhs_id]].get();
      auto head = nodes_[sequence_to_node[it.rhs_id] + 1 - it.strand].get();

      auto length = it.lhs_begin - it.rhs_begin;
      auto length_pair = (piles_[it.rhs_id]->length() - it.rhs_end) -
                         (piles_[it.lhs_id]->length() - it.lhs_end);

      if (it.score == 4) {
        std::swap(head, tail);
        length *= -1;
        length_pair *= -1;
      }

      auto edge = std::make_shared<Edge>(tail, head, length);
      edges_.emplace_back(edge);
      edges_.emplace_back(std::make_shared<Edge>(head->pair, tail->pair,
                                                 length_pair));  // NOLINT
      edge->pair = edges_.back().get();
      edge->pair->pair = edge.get();
    }

    std::cerr << "[raven::Graph::Construct] stored " << edges_.size()
              << " edges "  // NOLINT
              << std::fixed << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -4) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Construct] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  std::cerr << "[raven::Graph::Construct] " << std::fixed
            << timer.elapsed_time() << "s" << std::endl;
}  // NOLINT

std::size_t Graph::GreedyConstruct(
    std::vector<std::unique_ptr<biosoup::Sequence>>& sequences) {
  using detail::OverlapCategory;
  using detail::OverlapSide;

  sequences.resize(1);

  biosoup::Timer timer{};
  auto fillers = util::LoadFillerSeqs();

  auto const n_unitigs = sequences.size();
  auto n_fillers = fillers.size();

  auto const unitigs_begin = std::size_t{0};
  auto const unitigs_end = unitigs_begin + n_unitigs;

  auto const fillers_begin = unitigs_end;
  auto fillers_end = unitigs_end + n_fillers;

  std::unordered_set<std::uint32_t> relevant_fillers;
  std::unordered_map<std::uint32_t, std::uint32_t> node_indices;

  std::vector<detail::OverlapSideVecs> unitig_overlaps(n_unitigs);
  std::vector<std::vector<biosoup::Overlap>> overlaps;

  auto const is_filler_irrelevant =
      [&relevant_fillers](
          std::unique_ptr<biosoup::Sequence> const& filler) -> bool {
    return relevant_fillers.find(filler->id) == relevant_fillers.end();
  };

  auto const overlap_type =
      [&sequences](biosoup::Overlap const& ovlp) -> std::uint32_t {
    auto const& strand = ovlp.strand;

    std::uint32_t const lhs_len = sequences[ovlp.lhs_id]->data.size();
    std::uint32_t const& lhs_begin = ovlp.lhs_begin;
    std::uint32_t const& lhs_end = ovlp.lhs_end;

    std::uint32_t const rhs_len = sequences[ovlp.rhs_id]->data.size();
    std::uint32_t const rhs_begin =
        strand ? ovlp.rhs_begin : rhs_len - ovlp.rhs_end;
    std::uint32_t const rhs_end =
        strand ? ovlp.rhs_end : rhs_len - ovlp.rhs_begin;

    auto const overhang = std::min(lhs_begin, rhs_begin) +
                          std::min(lhs_len - lhs_end, rhs_len - rhs_end);

    if (lhs_end - lhs_begin < (lhs_end - lhs_begin + overhang) * 0.875 ||
        rhs_end - rhs_begin < (rhs_end - rhs_begin + overhang) * 0.875) {
      return 0;  // internal
    }
    if (lhs_begin <= rhs_begin && lhs_len - lhs_end <= rhs_len - rhs_end) {
      return 1;  // lhs contained
    }
    if (rhs_begin <= lhs_begin && rhs_len - rhs_end <= lhs_len - lhs_end) {
      return 2;  // rhs contained
    }
    if (lhs_begin > rhs_begin) {
      return 3;  // lhs -> rhs
    }
    return 4;  // rhs -> lhs
  };

  auto const overlap_set_score =
      [&overlap_type](biosoup::Overlap& ovlp) -> void {
    ovlp.score = overlap_type(ovlp);
  };

  // 0 -> irrelevant; 1 -> left, 2 -> right
  auto const overlap_category = [&](biosoup::Overlap& ovlp) -> OverlapCategory {
    overlap_set_score(ovlp);
    if (ovlp.score <= 2) {
      return OverlapCategory::kIrrelevant;
    }

    // unitigs always have a lower id
    auto const unitig_len =
        sequences[std::min(ovlp.lhs_id, ovlp.rhs_id)]->data.size();

    auto const left_delim = std::min(
        static_cast<std::size_t>(unitig_len * 0.005), std::size_t{5000});
    auto const right_delim =
        std::max(static_cast<std::size_t>(unitig_len * 0.995),
                 std::size_t{unitig_len - 5000});

    auto const unitig_ovlp_begin =
        ovlp.lhs_id < ovlp.rhs_id ? ovlp.lhs_begin : ovlp.rhs_begin;
    auto const unitig_ovlp_end =
        ovlp.lhs_id < ovlp.rhs_id ? ovlp.lhs_end : ovlp.rhs_end;

    auto const is_left_contained =
        unitig_ovlp_begin <= left_delim && unitig_ovlp_end <= left_delim;
    auto const is_right_contained =
        unitig_ovlp_end >= right_delim && unitig_ovlp_end >= right_delim;

    if (is_left_contained) {
      return OverlapCategory::kLeft;
    }

    if (is_right_contained) {
      return OverlapCategory::kRight;
    }

    return OverlapCategory::kIrrelevant;
  };

  auto const emplace_overlap = [&](biosoup::Overlap const& ovlp) -> void {
    if (ovlp.lhs_id < unitigs_end) {
      unitig_overlaps[ovlp.lhs_id].AddOverlap(sequences[ovlp.lhs_id], ovlp);
    } else if (ovlp.rhs_id < unitigs_end) {
      unitig_overlaps[ovlp.rhs_id].AddOverlap(sequences[ovlp.rhs_id],
                                              detail::OverlapReverse(ovlp));
    } else {
      overlaps[ovlp.lhs_id].emplace_back(ovlp);
      overlaps[ovlp.rhs_id].emplace_back(ovlp);
    }
  };

  auto const emplace_node = [&](std::unique_ptr<biosoup::Sequence> const& seq)
      -> std::shared_ptr<Node>& {
    nodes_.emplace_back(std::make_shared<Node>(*seq.get()));
    return nodes_.back();
  };

  auto const emplace_edge = [&](Node* tail, Node* head,
                                std::uint32_t len) -> std::shared_ptr<Edge>& {
    edges_.emplace_back(std::make_shared<Edge>(tail, head, len));
    return edges_.back();
  };

  // used in graph assembly construction
  auto const is_sequence_used = [&](std::uint32_t sequence_id) -> bool {
    return node_indices.find(sequence_id) != node_indices.end();
  };

  // returns a node index
  auto const node_from_sequence =
      [&](std::unique_ptr<biosoup::Sequence> const& sequence) -> std::uint32_t {
    auto const node_index = static_cast<std::uint32_t>(Node::num_objects);
    node_indices[sequence->id] = node_index;

    auto node_a = emplace_node(sequence);
    sequence->ReverseAndComplement();
    auto node_b = emplace_node(sequence);

    node_a->pair = node_b.get();
    node_b->pair = node_a.get();

    return node_index;
  };

  // returns an edge index
  auto const edge_from_overlap =
      [&](biosoup::Overlap const& ovlp) -> std::uint32_t {
    auto const edge_index = static_cast<std::uint32_t>(Edge::num_objects);

    auto const lhs_index = node_indices[ovlp.lhs_id];
    auto const rhs_index = node_indices[ovlp.rhs_id];

    auto tail = nodes_[lhs_index].get();
    auto head = nodes_[rhs_index].get();

    auto length = 1LL * ovlp.lhs_begin - ovlp.rhs_begin;
    auto length_pair = 1LL * ovlp.lhs_end - ovlp.rhs_end;

    if (ovlp.score == 4) {
      std::swap(tail, head);
      length_pair *= -1;
      length *= -1;
    }

    auto edge_a = emplace_edge(tail, head, length);
    auto edge_b = emplace_edge(head->pair, tail->pair, length_pair);

    edge_a->pair = edge_b.get();
    edge_b->pair = edge_a.get();

    return edge_index;
  };

  timer.Start();

  sequences = util::MergeSequences(sequences, fillers);
  util::NormalizeSeqIds(sequences);

  minimizer_engine_.Minimize(sequences.begin() + unitigs_begin,
                             sequences.begin() + unitigs_end);

  std::cerr << "[raven::Graph::GreedyConstruct] minimized " << n_unitigs
            << " unitigs " << timer.Stop() << 's' << std::endl;

  timer.Start();

  std::size_t overlap_batch_size = 0;
  overlaps.resize(sequences.size());
  std::vector<std::future<std::vector<biosoup::Overlap>>> ovlp_futures;
  for (std::size_t i = fillers_begin; i < fillers_end; ++i) {
    overlap_batch_size += sequences[i]->data.size();
    ovlp_futures.emplace_back(thread_pool_->Submit(
        [&](std::size_t const pos) -> std::vector<biosoup::Overlap> {
          return minimizer_engine_.Map(sequences[pos], true, false, true);
        },
        i));

    if (i + 1 != fillers_end && overlap_batch_size < constants::kOvlpBatchLim) {
      continue;
    }

    for (auto& ovlp_vec : ovlp_futures) {
      for (auto& ovlp : ovlp_vec.get()) {
        auto const category = overlap_category(ovlp);
        if (category != OverlapCategory::kIrrelevant) {
          relevant_fillers.emplace(ovlp.lhs_id);
        }
      }
    }

    ovlp_futures.clear();
  }

  sequences.erase(
      std::remove_if(sequences.begin() + fillers_begin,
                     sequences.begin() + fillers_end, is_filler_irrelevant),
      sequences.begin() + fillers_end);

  util::NormalizeSeqIds(sequences);

  fillers_end = sequences.size();
  n_fillers = fillers_end - fillers_begin;

  std::cerr << "[raven::GreedyConstruct] filtered out relevant fillers ("
            << n_fillers << ") " << timer.Stop() << "s" << std::endl;

  // map unitigs and valid fillers to earchother

  overlaps.resize(sequences.size());

  std::size_t sequence_batch_size = 0;
  for (std::size_t i = 0, j = 0; i < sequences.size(); ++i) {
    sequence_batch_size += sequences[i]->data.size();
    if (i + 1 != sequences.size() &&
        sequence_batch_size < constants::kSeqsBatchLim) {
      continue;
    }

    timer.Start();

    minimizer_engine_.Minimize(sequences.begin() + j,
                               sequences.begin() + i + 1);

    std::cerr << "[raven::Graph::GreedyConstruct] minimized " << j << " - "
              << i + 1 << " / " << sequences.size() << " " << std::fixed
              << timer.Stop() << "s" << std::endl;

    timer.Start();

    if (i >= fillers_begin) {
      std::size_t overlap_batch_size = 0;
      std::vector<std::future<std::vector<biosoup::Overlap>>> ovlp_futures;
      for (std::size_t k = j; k < i + 1; ++k) {
        overlap_batch_size += sequences[k]->data.size();
        ovlp_futures.emplace_back(thread_pool_->Submit(
            [&](std::size_t pos) -> std::vector<biosoup::Overlap> {
              return minimizer_engine_.Map(sequences[pos], true, true, true);
            },
            k));

        if (k != i && overlap_batch_size < constants::kOvlpBatchLim) {
          continue;
        }

        for (auto& ovlp_vec : ovlp_futures) {
          for (auto& ovlp : ovlp_vec.get()) {
            emplace_overlap(ovlp);
          }
        }
      }
    }

    std::cerr << "[raven::Graph::GreedyConstruct] mapped sequences "
              << std::fixed << timer.Stop() << "s" << std::endl;

    j = i + 1;
  }

  for (std::size_t i = unitigs_begin; i < unitigs_end; ++i) {
    overlaps[i] = unitig_overlaps[i].MergedSides();
  }

  timer.Start();

  // construct the graph from overlaps

  // for (auto const& sequence : sequences) {
  //   node_from_sequence(sequence);
  // }

  // for (auto const& ovlp_vec : overlaps) {
  //   for (auto const& ovlp : ovlp_vec) {
  //     edge_from_overlap(ovlp);
  //   }
  // }

  std::deque<std::uint32_t> ovlp_segments;
  auto const construction_step = [&](biosoup::Overlap const& ovlp) -> void {
    if (!is_sequence_used(ovlp.rhs_id)) {
      node_from_sequence(sequences[ovlp.rhs_id]);
      ovlp_segments.emplace_back(ovlp.rhs_id);
    }

    edge_from_overlap(ovlp);
  };

  for (std::size_t i = unitigs_begin; i < unitigs_end; ++i) {
    node_from_sequence(sequences[i]);
    auto const overlap_vec = overlaps[i];
    for (auto const& overlap : overlap_vec) {
      construction_step(overlap);
    }
  }

  while (!ovlp_segments.empty()) {
    auto const segment_id = ovlp_segments.front();
    ovlp_segments.pop_front();

    auto const overlap_vec = overlaps[segment_id];
    for (auto const& overlap : overlap_vec) {
      construction_step(overlap);
    }
  }

  std::cerr << "[raven::GreedyConstruct] assembly graf constructed "
            << timer.Stop() << 's' << std::endl;

  std::cerr << "[raven::GreedyConstruct] stored " << nodes_.size() << " nodes"
            << std::endl;

  std::cerr << "[raven::GreedyConstruct] stored " << edges_.size() << " edges"
            << std::endl;

  return n_unitigs;
}  // namespace raven

void Graph::Assemble() {
  if (stage_ < -3 || stage_ > -1) {
    return;
  }

  biosoup::Timer timer{};

  if (stage_ == -3) {  // remove transitive edges
    timer.Start();

    RemoveTransitiveEdges();

    std::cerr << "[raven::Graph::Assemble] removed transitive edges "
              << std::fixed << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -3) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -2) {  // remove tips and bubbles
    timer.Start();

    while (true) {
      std::uint32_t num_changes = RemoveTips();
      num_changes += RemoveBubbles();
      if (num_changes == 0) {
        break;
      }
    }

    std::cerr << "[raven::Graph::Assemble] removed tips and bubbles "
              << std::fixed << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -2) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -1) {  // remove long edges
    timer.Start();

    CreateUnitigs(42);  // speed up force directed layout
    RemoveLongEdges(16);

    std::cerr << "[raven::Graph::Assemble] removed long edges " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  if (stage_ == -1) {  // checkpoint
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Assemble] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }

  timer.Start();

  while (true) {  // TODO(rvaser): check if necessary
    std::uint32_t num_changes = RemoveTips();
    num_changes += RemoveBubbles();
    if (num_changes == 0) {
      break;
    }
  }

  timer.Stop();
  std::cerr << "[raven::Graph::Assemble] " << std::fixed << timer.elapsed_time()
            << "s" << std::endl;
}

void Graph::GreedyAssemble(std::size_t const n_expected) {
  std::unordered_set<std::uint32_t> valid_nodes;
  std::unordered_set<std::uint32_t> marked_edges;

  using detail::ExpandDir;

  biosoup::Timer timer;

  auto const greedy_expand = [&](Node const* starting_node,
                                 ExpandDir const dir) -> bool {
    std::unordered_set<std::uint32_t> dfs_visited;

    auto const not_visited = [&](Node const* node) -> bool {
      return valid_nodes.find(node->id) == valid_nodes.end() &&
             dfs_visited.find(node->id) == dfs_visited.end();
    };

    auto const mark_edge = [&](Edge const* edge_ptr) -> void {
      marked_edges.insert(edge_ptr->id);
    };

    auto const mark_edges_except = [&](std::vector<Edge*> edges,
                                       Edge const* excluded_edge) -> void {
      for (auto const edge_ptr : edges) {
        if (edge_ptr->id != excluded_edge->id) {
          mark_edge(edge_ptr);
          mark_edge(edge_ptr->pair);
        }
      }
    };

    using expand_fn_t = std::function<bool(Node const*)>;

    expand_fn_t const expand_left = [&](Node const* curr_node) -> bool {
      dfs_visited.insert(curr_node->id);
      for (auto const curr_in_edge : curr_node->inedges) {
        auto const nxt_node = curr_in_edge->tail;
        if (nxt_node->id == starting_node->id ||
            (not_visited(nxt_node) && expand_left(nxt_node))) {
          valid_nodes.insert(nxt_node->id);
          mark_edges_except(curr_node->inedges, curr_in_edge);
          return true;
        }
      }

      return false;
    };

    expand_fn_t const expand_right = [&](Node const* curr_node) -> bool {
      dfs_visited.insert(curr_node->id);
      for (auto const curr_out_edge : curr_node->outedges) {
        auto const nxt_node = curr_out_edge->head;
        if (nxt_node->id == starting_node->id ||
            (not_visited(nxt_node) && expand_right(nxt_node))) {
          valid_nodes.insert(nxt_node->id);
          mark_edges_except(curr_node->outedges, curr_out_edge);

          return true;
        }
      }

      return false;
    };

    expand_fn_t const expand_fn = [&]() -> expand_fn_t const {
      if (dir == ExpandDir::kLeft) {
        return expand_left;
      } else {
        return expand_right;
      }
    }();

    if (expand_fn(starting_node)) {
      valid_nodes.insert(starting_node->id);
      return true;
    }

    return false;
  };

  auto const sort_edges_by_len = [](std::vector<Edge*> edges) -> void {
    std::sort(edges.begin(), edges.end(),
              [](Edge* const a, Edge* const b) -> bool {
                return a->length > b->length;
              });
  };

  auto const sort_node_edges_by_len = [&sort_edges_by_len](Node* node) -> void {
    sort_edges_by_len(node->inedges);
    sort_edges_by_len(node->outedges);
  };

  timer.Start();
  for (auto& node : nodes_) {
    sort_node_edges_by_len(node.get());
  }

  std::cerr << "[raven::Graph::GreedyAssemble] sorted edges by len "
            << timer.Stop() << 's' << std::endl;

  timer.Start();
  for (std::size_t i = 0; i < n_expected * 2; i += 2) {
    auto const& curr_node = nodes_[i];
    std::cerr << "[raven::de] starting from: " << curr_node->name << std::endl;
    if (greedy_expand(curr_node.get(), ExpandDir::kLeft) ||
        greedy_expand(curr_node.get(), ExpandDir::kRight)) {
      std::cerr << "[raven::de] found path from: " << curr_node->name
                << std::endl;
      RemoveEdges(marked_edges);
      marked_edges.clear();
    }
  }

  std::cerr << "[raven::Graph::GreedyAssemble] found paths " << timer.Stop()
            << 's' << std::endl;

  std::cerr << "[raven::Graph::GreedyAssemble] " << timer.elapsed_time() << 's'
            << std::endl;
}

std::uint32_t Graph::RemoveTransitiveEdges() {
  auto is_comparable = [](double const a, double const b) -> bool {
    double constexpr eps = 0.12;
    return (a >= b * (1 - eps) && a <= b * (1 + eps)) ||
           (b >= a * (1 - eps) && b <= a * (1 + eps));
  };

  std::vector<Edge*> candidate(nodes_.size(), nullptr);
  std::unordered_set<std::uint32_t> marked_edges;
  for (const auto& it : nodes_) {
    if (it == nullptr) {
      continue;
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = jt;
    }
    for (auto jt : it->outedges) {
      for (auto kt : jt->head->outedges) {
        if (candidate[kt->head->id] &&
            is_comparable(jt->length + kt->length,
                          candidate[kt->head->id]->length)) {  // NOLINT
          marked_edges.emplace(candidate[kt->head->id]->id);
          marked_edges.emplace(candidate[kt->head->id]->pair->id);
        }
      }
    }
    for (auto jt : it->outedges) {
      candidate[jt->head->id] = nullptr;
    }
  }

  for (auto i : marked_edges) {  // store for force directed layout
    if (i & 1) {
      auto lhs = edges_[i]->tail->id & ~1UL;
      auto rhs = edges_[i]->head->id & ~1UL;
      nodes_[lhs]->transitive.emplace(rhs);
      nodes_[rhs]->transitive.emplace(lhs);
    }
  }

  RemoveEdges(marked_edges);
  return marked_edges.size() / 2;
}

std::uint32_t Graph::RemoveTips() {
  std::uint32_t num_tips = 0;
  std::vector<char> is_visited(nodes_.size(), 0);

  for (const auto& it : nodes_) {
    if (it == nullptr || is_visited[it->id] || !it->is_tip()) {
      continue;
    }
    bool is_circular = false;
    std::uint32_t num_sequences = 0;

    auto end = it.get();
    while (!end->is_junction()) {
      num_sequences += end->count;
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 || end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    if (is_circular || end->outdegree() == 0 || num_sequences > 5) {
      continue;
    }

    std::unordered_set<std::uint32_t> marked_edges;
    for (auto jt : end->outedges) {
      if (jt->head->indegree() > 1) {
        marked_edges.emplace(jt->id);
        marked_edges.emplace(jt->pair->id);
      }
    }
    if (marked_edges.size() / 2 == end->outedges.size()) {  // delete whole
      auto begin = it.get();
      while (begin != end) {
        marked_edges.emplace(begin->outedges.front()->id);
        marked_edges.emplace(begin->outedges.front()->pair->id);
        begin = begin->outedges.front()->head;
      }
      ++num_tips;
    }
    RemoveEdges(marked_edges, true);
  }

  return num_tips;
}

std::uint32_t Graph::RemoveBubbles() {
  std::vector<std::uint32_t> distance(nodes_.size(), 0);
  std::vector<std::uint32_t> n_nodes_to(nodes_.size(), 0);  // TODO: Remove?
  std::vector<Node*> predecessor(nodes_.size(), nullptr);

  // path helper functions
  auto path_extract = [&](Node* begin, Node* end) -> std::vector<Node*> {
    std::vector<Node*> dst;
    while (end != begin) {
      dst.emplace_back(end);
      end = predecessor[end->id];
    }
    dst.emplace_back(begin);
    std::reverse(dst.begin(), dst.end());
    return dst;
  };

  auto path_type = [](const std::vector<Node*>& path) -> bool {
    if (path.empty()) {
      return false;
    }
    for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
      if (path[i]->is_junction()) {
        return false;  // complex
      }
    }
    return true;  // without branches
  };

  auto bubble_type = [&](const std::vector<Node*>& lhs,
                         const std::vector<Node*>& rhs) -> bool {
    if (lhs.empty() || rhs.empty()) {
      return false;
    }
    std::unordered_set<Node*> intersection;
    for (auto it : lhs) {
      intersection.emplace(it);
    }
    for (auto it : rhs) {
      intersection.emplace(it);
    }
    if (lhs.size() + rhs.size() - 2 != intersection.size()) {
      return false;
    }
    for (auto it : lhs) {
      if (intersection.count(it->pair) != 0) {
        return false;
      }
    }
    if (path_type(lhs) && path_type(rhs)) {  // both without branches
      return true;
    }

    auto path_sequence = [](const std::vector<Node*>& path)
        -> std::unique_ptr<biosoup::Sequence> {
      auto sequence =
          std::unique_ptr<biosoup::Sequence>(new biosoup::Sequence());
      for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
        for (auto it : path[i]->outedges) {
          if (it->head == path[i + 1]) {
            sequence->data += it->Label();
            break;
          }
        }
      }
      sequence->data += path.back()->data;
      return std::move(sequence);
    };

    auto ls = path_sequence(lhs);
    auto rs = path_sequence(rhs);

    if (std::min(ls->data.size(), rs->data.size()) <
        std::max(ls->data.size(), rs->data.size()) * 0.8) {
      return false;
    }

    auto overlaps = minimizer_engine_.Map(ls, rs);

    std::uint32_t matches = 0;
    for (const auto& it : overlaps) {
      matches = std::max(matches, it.score);
    }
    return matches > 0.5 * std::min(ls->data.size(), rs->data.size());
  };
  // path helper functions

  std::uint32_t num_bubbles = 0;
  for (const auto& it : nodes_) {
    if (it == nullptr || it->outdegree() < 2) {
      continue;
    }

    // BFS
    Node* begin = it.get();
    Node* end = nullptr;
    Node* other_end = nullptr;
    std::deque<Node*> que{begin};
    std::vector<Node*> visited{1, begin};
    while (!que.empty() && !end) {
      auto jt = que.front();
      que.pop_front();

      for (auto kt : jt->outedges) {
        if (kt->head == begin) {  // cycle
          continue;
        }
        // if (distance[jt->id] + kt->length > 500000) {  // out of reach
        //   continue;
        // }
        if (n_nodes_to[jt->id] > 3400) {  // out of reach
          continue;
        }
        // distance[kt->head->id] = distance[jt->id] + kt->length;
        n_nodes_to[kt->head->id] = n_nodes_to[jt->id] + 1;
        visited.emplace_back(kt->head);
        que.emplace_back(kt->head);

        if (predecessor[kt->head->id]) {  // found bubble
          end = kt->head;
          other_end = jt;
          break;
        }

        predecessor[kt->head->id] = jt;
      }
    }

    std::unordered_set<std::uint32_t> marked_edges;
    if (end) {
      auto lhs = path_extract(begin, end);
      auto rhs = path_extract(begin, other_end);
      rhs.emplace_back(end);

      if (bubble_type(lhs, rhs)) {
        std::uint32_t lhs_count = 0;
        for (auto jt : lhs) {
          lhs_count += jt->count;
        }
        std::uint32_t rhs_count = 0;
        for (auto jt : rhs) {
          rhs_count += jt->count;
        }
        marked_edges = FindRemovableEdges(lhs_count > rhs_count ? rhs : lhs);
        if (marked_edges.empty()) {
          marked_edges = FindRemovableEdges(lhs_count > rhs_count ? lhs : rhs);
        }
      }
    }

    for (auto jt : visited) {
      // distance[jt->id] = 0;
      n_nodes_to[jt->id] = 0;
      predecessor[jt->id] = nullptr;
    }

    RemoveEdges(marked_edges, true);
    num_bubbles += 1 - marked_edges.empty();
  }

  return num_bubbles;
}

std::uint32_t Graph::RemoveLongEdges(std::uint32_t num_rounds) {
  std::uint32_t num_long_edges = 0;

  for (std::uint32_t i = 0; i < num_rounds; ++i) {
    CreateForceDirectedLayout();

    std::unordered_set<std::uint32_t> marked_edges;
    for (const auto& it : nodes_) {
      if (it == nullptr || it->outdegree() < 2) {
        continue;
      }
      for (auto jt : it->outedges) {
        for (auto kt : it->outedges) {
          if (jt != kt && jt->weight * 2.0 < kt->weight) {  // TODO(rvaser)
            marked_edges.emplace(kt->id);
            marked_edges.emplace(kt->pair->id);
          }
        }
      }
    }
    RemoveEdges(marked_edges);
    num_long_edges += marked_edges.size() / 2;

    RemoveTips();
  }

  return num_long_edges;
}

void Graph::CreateForceDirectedLayout(const std::string& path) {
  std::ofstream os;
  bool is_first = true;
  if (!path.empty()) {
    os.open(path);
    os << "{" << std::endl;
  }

  std::vector<std::unordered_set<std::uint32_t>> components;
  std::vector<char> is_visited(nodes_.size(), 0);
  for (std::uint32_t i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i] == nullptr || is_visited[i]) {
      continue;
    }

    components.resize(components.size() + 1);

    std::deque<std::uint32_t> que = {i};
    while (!que.empty()) {
      std::uint32_t j = que.front();
      que.pop_front();

      if (is_visited[j]) {
        continue;
      }
      const auto& node = nodes_[j];
      is_visited[node->id] = 1;
      is_visited[node->pair->id] = 1;
      components.back().emplace((node->id >> 1) << 1);

      for (auto it : node->inedges) {
        que.emplace_back(it->tail->id);
      }
      for (auto it : node->outedges) {
        que.emplace_back(it->head->id);
      }
    }
  }
  std::vector<char>().swap(is_visited);

  std::sort(components.begin(), components.end(),
            [](const std::unordered_set<std::uint32_t>& lhs,
               const std::unordered_set<std::uint32_t>& rhs) {
              return lhs.size() > rhs.size();
            });

  static std::uint64_t seed = 21;
  seed <<= 1;

  std::mt19937 generator(seed);
  std::uniform_real_distribution<> distribution(0., 1.);

  struct Point {
    Point() = default;
    Point(double x, double y) : x(x), y(y) {}

    bool operator==(const Point& other) const {
      return x == other.x && y == other.y;
    }
    Point operator+(const Point& other) const {
      return Point(x + other.x, y + other.y);
    }
    Point& operator+=(const Point& other) {
      x += other.x;
      y += other.y;
      return *this;
    }
    Point operator-(const Point& other) const {
      return Point(x - other.x, y - other.y);
    }
    Point operator*(double c) const { return Point(x * c, y * c); }
    Point& operator/=(double c) {
      x /= c;
      y /= c;
      return *this;
    }
    double norm() const { return sqrt(x * x + y * y); }

    double x;
    double y;
  };

  struct Quadtree {
    Quadtree(Point nucleus, double width)
        : nucleus(nucleus), width(width), center(0, 0), mass(0), subtrees() {}

    bool add(const Point& p) {
      if (nucleus.x - width > p.x || p.x > nucleus.x + width ||
          nucleus.y - width > p.y || p.y > nucleus.y + width) {
        return false;
      }
      ++mass;
      if (mass == 1) {
        center = p;
      } else if (subtrees.empty()) {
        if (center == p) {
          return true;
        }
        double w = width / 2;
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y + w), w);
        subtrees.emplace_back(Point(nucleus.x - w, nucleus.y - w), w);
        subtrees.emplace_back(Point(nucleus.x + w, nucleus.y - w), w);
        for (auto& it : subtrees) {
          if (it.add(center)) {
            break;
          }
        }
      }
      for (auto& it : subtrees) {
        if (it.add(p)) {
          break;
        }
      }
      return true;
    }

    void centre() {
      if (subtrees.empty()) {
        return;
      }
      center = Point(0, 0);
      for (auto& it : subtrees) {
        it.centre();
        center += it.center * it.mass;
      }
      center /= mass;
    }

    Point force(const Point& p, double k) const {
      auto delta = p - center;
      auto distance = delta.norm();
      if (width * 2 / distance < 1) {
        return delta * (mass * (k * k) / (distance * distance));
      }
      delta = Point(0, 0);
      for (const auto& it : subtrees) {
        delta += it.force(p, k);
      }
      return delta;
    }

    Point nucleus;
    double width;
    Point center;
    std::uint32_t mass;
    std::vector<Quadtree> subtrees;
  };

  std::uint32_t c = 0;
  for (const auto& component : components) {
    if (component.size() < 6) {
      continue;
    }

    bool has_junctions = false;
    for (const auto& it : component) {
      if (nodes_[it]->is_junction()) {
        has_junctions = true;
        break;
      }
    }
    if (has_junctions == false) {
      continue;
    }

    // update transitive edges
    for (const auto& n : component) {
      std::unordered_set<std::uint32_t> valid;
      for (const auto& m : nodes_[n]->transitive) {
        if (component.find(m) != component.end()) {
          valid.emplace(m);
        }
      }
      nodes_[n]->transitive.swap(valid);
    }

    std::uint32_t num_iterations = 100;
    double k = sqrt(1. / static_cast<double>(component.size()));
    double t = 0.1;
    double dt = t / static_cast<double>(num_iterations + 1);

    std::vector<Point> points(nodes_.size());
    for (const auto& it : component) {
      points[it].x = distribution(generator);
      points[it].y = distribution(generator);
    }

    for (std::uint32_t i = 0; i < num_iterations; ++i) {
      Point x = {0, 0}, y = {0, 0};
      for (const auto& n : component) {
        x.x = std::min(x.x, points[n].x);
        x.y = std::max(x.y, points[n].x);
        y.x = std::min(y.x, points[n].y);
        y.y = std::max(y.y, points[n].y);
      }
      double w = (x.y - x.x) / 2, h = (y.y - y.x) / 2;

      Quadtree tree(Point(x.x + w, y.x + h), std::max(w, h) + 0.01);
      for (const auto& n : component) {
        tree.add(points[n]);
      }
      tree.centre();

      std::vector<std::future<void>> thread_futures;
      std::vector<Point> displacements(nodes_.size(), Point(0, 0));

      auto thread_task = [&](std::uint32_t n) -> void {
        auto displacement = tree.force(points[n], k);
        for (auto e : nodes_[n]->inedges) {
          auto m = (e->tail->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (auto e : nodes_[n]->outedges) {
          auto m = (e->head->id >> 1) << 1;
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        for (const auto& m : nodes_[n]->transitive) {
          auto delta = points[n] - points[m];
          auto distance = delta.norm();
          if (distance < 0.01) {
            distance = 0.01;
          }
          displacement += delta * (-1. * distance / k);
        }
        auto length = displacement.norm();
        if (length < 0.01) {
          length = 0.1;
        }
        displacements[n] = displacement * (t / length);
        return;
      };

      for (const auto& n : component) {
        thread_futures.emplace_back(thread_pool_->Submit(thread_task, n));
      }
      for (const auto& it : thread_futures) {
        it.wait();
      }
      for (const auto& n : component) {
        points[n] += displacements[n];
      }

      t -= dt;
    }

    for (const auto& it : edges_) {
      if (it == nullptr || it->id & 1) {
        continue;
      }
      auto n = (it->tail->id >> 1) << 1;
      auto m = (it->head->id >> 1) << 1;

      if (component.find(n) != component.end() &&
          component.find(m) != component.end()) {
        it->weight = (points[n] - points[m]).norm();
        it->pair->weight = it->weight;
      }
    }

    if (!path.empty()) {
      if (!is_first) {
        os << "," << std::endl;
      }
      is_first = false;

      os << "    \"component_" << c++ << "\": {" << std::endl;

      bool is_first_node = true;
      os << "      \"nodes\": {" << std::endl;
      for (const auto& it : component) {
        if (!is_first_node) {
          os << "," << std::endl;
        }
        is_first_node = false;
        os << "        \"" << it << "\": [";
        os << points[it].x << ", ";
        os << points[it].y << ", ";
        os << (nodes_[it]->is_junction() ? 1 : 0) << ", ";
        os << nodes_[it]->count << "]";
      }
      os << std::endl << "      }," << std::endl;

      bool is_first_edge = true;
      os << "      \"edges\": [" << std::endl;
      for (const auto& it : component) {
        for (auto e : nodes_[it]->inedges) {
          auto o = (e->tail->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (auto e : nodes_[it]->outedges) {
          auto o = (e->head->id >> 1) << 1;
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 0]";
        }
        for (const auto& o : nodes_[it]->transitive) {
          if (it < o) {
            continue;
          }
          if (!is_first_edge) {
            os << "," << std::endl;
          }
          is_first_edge = false;
          os << "        [\"" << it << "\", \"" << o << "\", 1]";
        }
      }
      os << std::endl << "      ]" << std::endl;
      os << "    }";
    }
  }

  if (!path.empty()) {
    os << std::endl << "}";
    os << std::endl;
    os.close();
  }
}

void Graph::Polish(
    const std::vector<std::unique_ptr<biosoup::Sequence>>& sequences,
    std::uint8_t match, std::uint8_t mismatch, std::uint8_t gap,
    std::uint32_t cuda_poa_batches, bool cuda_banded_alignment,
    std::uint32_t cuda_alignment_batches, std::uint32_t num_rounds) {
  if (sequences.empty() || num_rounds == 0) {
    return;
  }

  auto unitigs = GetUnitigs();
  if (unitigs.empty()) {
    return;
  }

  double q = 0.;
  for (const auto& it : sequences) {
    if (it->quality.empty()) {
      continue;
    }
    double p = 0.;
    for (const auto& jt : it->quality) {
      p += jt - 33;
    }
    q += p / it->quality.size();
  }
  if (q == 0.) {  // when all values equal to '!'
    for (const auto& it : sequences) {
      it->quality.clear();
    }
  } else {
    q /= sequences.size();
  }

  auto polisher = racon::Polisher::Create(
      q, 0.3, 500, true, match, mismatch, gap, thread_pool_, cuda_poa_batches,
      cuda_banded_alignment, cuda_alignment_batches);

  while (stage_ < static_cast<std::int32_t>(num_rounds)) {
    polisher->Initialize(unitigs, sequences);
    auto polished = polisher->Polish(false);
    unitigs.swap(polished);

    for (const auto& it : unitigs) {  // store unitigs
      const auto& node = nodes_[std::atoi(&it->name[3])];
      std::size_t tag;
      if ((tag = it->name.rfind(':')) != std::string::npos) {
        if (std::atof(&it->name[tag + 1]) > 0) {
          node->is_polished = true;
          node->data = it->data;
          it->ReverseAndComplement();
          node->pair->data = it->data;
        }
      }
    }

    biosoup::Timer timer{};
    timer.Start();

    ++stage_;
    Store();

    std::cerr << "[raven::Graph::Polish] reached checkpoint " << std::fixed
              << timer.Stop() << "s" << std::endl;
  }
}

std::uint32_t Graph::CreateUnitigs(std::uint32_t epsilon) {
  std::unordered_set<std::uint32_t> marked_edges;
  std::vector<std::shared_ptr<Node>> unitigs;
  std::vector<std::shared_ptr<Edge>> unitig_edges;
  std::vector<std::uint32_t> node_updates(nodes_.size(), 0);
  std::vector<char> is_visited(nodes_.size(), 0);

  for (const auto& it : nodes_) {
    if (it == nullptr || is_visited[it->id] || it->is_junction()) {
      continue;
    }

    std::uint32_t extension = 1;

    bool is_circular = false;
    auto begin = it.get();
    while (!begin->is_junction()) {  // extend left
      is_visited[begin->id] = 1;
      is_visited[begin->pair->id] = 1;
      if (begin->indegree() == 0 ||
          begin->inedges.front()->tail->is_junction()) {
        break;
      }
      begin = begin->inedges.front()->tail;
      ++extension;
      if (begin == it.get()) {
        is_circular = true;
        break;
      }
    }

    auto end = it.get();
    while (!end->is_junction()) {  // extend right
      is_visited[end->id] = 1;
      is_visited[end->pair->id] = 1;
      if (end->outdegree() == 0 || end->outedges.front()->head->is_junction()) {
        break;
      }
      end = end->outedges.front()->head;
      ++extension;
      if (end == it.get()) {
        is_circular = true;
        break;
      }
    }

    // case: node is a junction
    if (!is_circular && begin == end) {
      continue;
    }

    if (!is_circular && extension < 2 * epsilon + 2) {
      continue;
    }

    if (begin != end) {
      for (std::uint32_t i = 0; i < epsilon;
           ++i) {  // remove nodes near junctions
        begin = begin->outedges.front()->head;
      }
      for (std::uint32_t i = 0; i < epsilon; ++i) {
        end = end->inedges.front()->tail;
      }
    }

    auto unitig = std::make_shared<Node>(begin, end);
    unitigs.emplace_back(unitig);
    unitigs.emplace_back(std::make_shared<Node>(end->pair, begin->pair));
    unitig->pair = unitigs.back().get();
    unitig->pair->pair = unitig.get();

    if (begin != end) {  // connect unitig to graph
      if (begin->indegree()) {
        marked_edges.emplace(begin->inedges.front()->id);
        marked_edges.emplace(begin->inedges.front()->pair->id);

        auto edge =
            std::make_shared<Edge>(begin->inedges.front()->tail, unitig.get(),
                                   begin->inedges.front()->length);
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Edge>(
            unitig->pair, begin->inedges.front()->pair->head,
            begin->inedges.front()->pair->length + unitig->pair->data.size() -
                begin->pair->data.size()));  // NOLINT
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
      if (end->outdegree()) {
        marked_edges.emplace(end->outedges.front()->id);
        marked_edges.emplace(end->outedges.front()->pair->id);

        auto edge = std::make_shared<Edge>(
            unitig.get(), end->outedges.front()->head,
            end->outedges.front()->length + unitig->data.size() -
                end->data.size());  // NOLINT
        unitig_edges.emplace_back(edge);
        unitig_edges.emplace_back(std::make_shared<Edge>(
            end->outedges.front()->pair->tail, unitig->pair,
            end->outedges.front()->pair->length));
        edge->pair = unitig_edges.back().get();
        edge->pair->pair = edge.get();
      }
    }

    auto jt = begin;
    while (true) {
      marked_edges.emplace(jt->outedges.front()->id);
      marked_edges.emplace(jt->outedges.front()->pair->id);

      // update transitive edges
      node_updates[jt->id & ~1UL] = unitig->id;
      unitig->transitive.insert(nodes_[jt->id & ~1UL]->transitive.begin(),
                                nodes_[jt->id & ~1UL]->transitive.end());

      if ((jt = jt->outedges.front()->head) == end) {
        break;
      }
    }
  }

  nodes_.insert(nodes_.end(), unitigs.begin(), unitigs.end());
  edges_.insert(edges_.end(), unitig_edges.begin(), unitig_edges.end());
  RemoveEdges(marked_edges, true);

  for (const auto& it : nodes_) {  // update transitive edges
    if (it) {
      std::unordered_set<std::uint32_t> valid;
      for (auto jt : it->transitive) {
        valid.emplace(node_updates[jt] == 0 ? jt : node_updates[jt]);
      }
      it->transitive.swap(valid);
    }
  }

  return unitigs.size() / 2;
}

std::vector<std::unique_ptr<biosoup::Sequence>> Graph::GetUnitigs(
    bool drop_unpolished) {
  CreateUnitigs();

  biosoup::Sequence::num_objects = 0;

  std::vector<std::unique_ptr<biosoup::Sequence>> dst;
  for (const auto& it : nodes_) {
    if (it == nullptr || it->is_rc() || !it->is_unitig()) {
      continue;
    }
    if (drop_unpolished && !it->is_polished) {
      continue;
    }

    std::string name = it->name + " LN:i:" + std::to_string(it->data.size()) +
                       " RC:i:" + std::to_string(it->count) +
                       " XO:i:" + std::to_string(it->is_circular);

    dst.emplace_back(new biosoup::Sequence(name, it->data));
  }

  return dst;
}

void Graph::RemoveEdges(const std::unordered_set<std::uint32_t>& indices,
                        bool remove_nodes) {
  auto erase_remove = [](std::vector<Edge*>& edges, Edge* marked) -> void {
    edges.erase(std::remove(edges.begin(), edges.end(), marked), edges.end());
  };

  std::unordered_set<std::uint32_t> node_indices;
  for (auto i : indices) {
    if (remove_nodes) {
      node_indices.emplace(edges_[i]->tail->id);
      node_indices.emplace(edges_[i]->head->id);
    }
    erase_remove(edges_[i]->tail->outedges, edges_[i].get());
    erase_remove(edges_[i]->head->inedges, edges_[i].get());
  }
  if (remove_nodes) {
    for (auto i : node_indices) {
      if (nodes_[i]->outdegree() == 0 && nodes_[i]->indegree() == 0) {
        nodes_[i].reset();
      }
    }
  }
  for (auto i : indices) {
    edges_[i].reset();
  }
}

std::unordered_set<std::uint32_t> Graph::FindRemovableEdges(
    const std::vector<Node*>& path) {
  if (path.empty()) {
    return std::unordered_set<std::uint32_t>{};
  }

  auto find_edge = [](Node* tail, Node* head) -> Edge* {
    for (auto it : tail->outedges) {
      if (it->head == head) {
        return it;
      }
    }
    return nullptr;
  };

  // find first node with multiple in edges
  std::int32_t pref = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->indegree() > 1) {
      pref = i;
      break;
    }
  }

  // find last node with multiple out edges
  std::int32_t suff = -1;
  for (std::uint32_t i = 1; i < path.size() - 1; ++i) {
    if (path[i]->outdegree() > 1) {
      suff = i;
    }
  }

  std::unordered_set<std::uint32_t> dst;
  if (pref == -1 && suff == -1) {  // remove whole path
    for (std::uint32_t i = 0; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
    return dst;
  }

  if (pref != -1 && path[pref]->outdegree() > 1) {  // complex path
    return dst;                                     // empty
  }
  if (suff != -1 && path[suff]->indegree() > 1) {  // complex path
    return dst;                                    // empty
  }

  if (pref == -1) {  // remove everything after last suffix node
    for (std::uint32_t i = suff; i < path.size() - 1; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff == -1) {  // remove everything before first prefix node
    for (std::int32_t i = 0; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  } else if (suff < pref) {  // remove everything in between
    for (std::int32_t i = suff; i < pref; ++i) {
      auto it = find_edge(path[i], path[i + 1]);
      dst.emplace(it->id);
      dst.emplace(it->pair->id);
    }
  }
  return dst;  // empty
}

void Graph::PrintJSON(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  cereal::JSONOutputArchive archive(os);
  for (const auto& it : piles_) {
    if (it->is_invalid()) {
      continue;
    }
    archive(cereal::make_nvp(std::to_string(it->id()), *(it.get())));
  }
}

void Graph::PrintCSV(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size() << " RC:i:" << it->count << ","
       << it->pair->id << " [" << it->pair->id / 2 << "]"
       << " LN:i:" << it->pair->data.size() << " RC:i:" << it->pair->count
       << ",0,-" << std::endl;
  }
  for (const auto& it : edges_) {
    if (it == nullptr) {
      continue;
    }
    os << it->tail->id << " [" << it->tail->id / 2 << "]"
       << " LN:i:" << it->tail->data.size() << " RC:i:" << it->tail->count
       << "," << it->head->id << " [" << it->head->id / 2 << "]"
       << " LN:i:" << it->head->data.size() << " RC:i:" << it->head->count
       << ",1," << it->id << " " << it->length << " " << it->weight
       << std::endl;
  }
  for (const auto& it : nodes_) {  // circular edges TODO(rvaser): check
    if (it == nullptr || !it->is_circular) {
      continue;
    }
    os << it->id << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size() << " RC:i:" << it->count << "," << it->id
       << " [" << it->id / 2 << "]"
       << " LN:i:" << it->data.size() << " RC:i:" << it->count << ",1,-"
       << std::endl;
  }
  os.close();
}

void Graph::PrintGFA(const std::string& path) const {
  if (path.empty()) {
    return;
  }

  std::ofstream os(path);
  for (const auto& it : nodes_) {
    if ((it == nullptr) || it->is_rc() ||
        (it->count == 1 && it->outdegree() == 0 && it->indegree() == 0)) {
      continue;
    }
    os << "S\t" << it->name << "\t" << it->data << "\tLN:i:" << it->data.size()
       << "\tRC:i:" << it->count << std::endl;
    if (it->is_circular) {
      os << "L\t" << it->name << "\t" << '+' << "\t" << it->name << "\t" << '+'
         << "\t0M" << std::endl;
    }
  }
  for (const auto& it : edges_) {
    if (it == nullptr) {
      continue;
    }
    os << "L\t" << it->tail->name << "\t" << (it->tail->is_rc() ? '-' : '+')
       << "\t" << it->head->name << "\t" << (it->head->is_rc() ? '-' : '+')
       << "\t" << it->tail->data.size() - it->length << 'M' << std::endl;
  }
  os.close();
}

void Graph::Store() const {
  std::ofstream os("raven.cereal");
  try {
    cereal::BinaryOutputArchive archive(os);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Store] error: unable to store archive");
  }
}

void Graph::Load() {
  std::ifstream is("raven.cereal");
  try {
    cereal::BinaryInputArchive archive(is);
    archive(*this);
  } catch (std::exception&) {
    throw std::logic_error(
        "[raven::Graph::Load] error: unable to load archive");
  }
}

void Graph::Clear() {
  piles_.clear();
  nodes_.clear();
  edges_.clear();

  stage_ = -5;

  // reset global counters
  biosoup::Sequence::num_objects = 0;
  Node::num_objects = 0;
  Edge::num_objects = 0;
}

}  // namespace raven
