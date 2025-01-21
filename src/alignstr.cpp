/*
* The MIT License (MIT)
 *
 * Copyright (c) 2025 Yaroslav Riabtsev <yaroslav.riabtsev@rwth-aachen.de>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "alignstr.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

alignment::alignment(std::string seq_a, std::string seq_b, const score_matrix sm, const int d) : d(d), sm(sm),
    seq1(std::move(seq_a)),
    seq2(std::move(seq_b)),
    n(seq1.size()),
    m(seq2.size()) {
}

void alignment::needleman_wunsch() {
    if (!needleman_wunsch_matrix.empty()) {
        throw std::runtime_error("Needleman-Wunsch matrix is already computed");
    }
    needleman_wunsch_matrix.resize(m + 1, std::vector(n + 1, 0));
    std::vector sources(m + 1, std::vector(n + 1, 0));
    for (int i = 1; i <= m; ++i) {
        needleman_wunsch_matrix[i][0] = -i * d;
        sources[i][0] = static_cast<int>(align_sources::delete_align);
    }
    for (int j = 1; j <= n; ++j) {
        needleman_wunsch_matrix[0][j] = -j * d;
        sources[0][j] = static_cast<int>(align_sources::insert_align);
    }
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            needleman_wunsch_matrix[i][j] = needleman_wunsch_matrix[i - 1][j] - d;
            sources[i][j] = static_cast<int>(align_sources::delete_align);

            if (const int _insert = needleman_wunsch_matrix[i][j - 1] - d; _insert > needleman_wunsch_matrix[i][j]) {
                needleman_wunsch_matrix[i][j] = _insert;
                sources[i][j] = static_cast<int>(align_sources::insert_align);
            } else if (_insert == needleman_wunsch_matrix[i][j]) {
                sources[i][j] |= static_cast<int>(align_sources::insert_align);
            }

            if (const int _match = needleman_wunsch_matrix[i - 1][j - 1] + similarity(seq1[j - 1], seq2[i - 1]);
                _match > needleman_wunsch_matrix[i][j]) {
                needleman_wunsch_matrix[i][j] = _match;
                sources[i][j] = static_cast<int>(align_sources::match_align);
            } else if (_match == needleman_wunsch_matrix[i][j]) {
                sources[i][j] |= static_cast<int>(align_sources::match_align);
            }
        }
    }

    // todo: backtrace
    // tip: needleman_wunsch() -> needleman_wunsch(std::vector<std::pair<std::string, std::string>>& best_align_pairs)
}

void alignment::print_needleman_wunsch() {
    if (needleman_wunsch_matrix.empty()) {
        needleman_wunsch();
    }
    std::cout << "Needleman-Wunsch matrix:" << std::endl;
    print(needleman_wunsch_matrix);
}

void alignment::smith_waterman() {
    if (!smith_waterman_matrix.empty()) {
        throw std::runtime_error("Smith-Waterman matrix is already computed");
    }
    smith_waterman_matrix.resize(m + 1, std::vector(n + 1, 0));
    std::vector sources(m + 1, std::vector(n + 1, 0));
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (const int _delete = smith_waterman_matrix[i - 1][j] - d;
                _delete > smith_waterman_matrix[i][j]) {
                smith_waterman_matrix[i][j] = _delete;
                sources[i][j] = static_cast<int>(align_sources::delete_align);
            } else if (_delete == smith_waterman_matrix[i][j]) {
                sources[i][j] |= static_cast<int>(align_sources::delete_align);
            }

            if (const int _insert = smith_waterman_matrix[i][j - 1] - d;
                _insert > smith_waterman_matrix[i][j]) {
                smith_waterman_matrix[i][j] = _insert;
                sources[i][j] = static_cast<int>(align_sources::insert_align);
            } else if (_insert == smith_waterman_matrix[i][j]) {
                sources[i][j] |= static_cast<int>(align_sources::insert_align);
            }

            if (const int _match = smith_waterman_matrix[i - 1][j - 1] + similarity(seq1[j - 1], seq2[i - 1]);
                _match > smith_waterman_matrix[i][j]) {
                smith_waterman_matrix[i][j] = _match;
                sources[i][j] = static_cast<int>(align_sources::match_align);
            } else if (_match == smith_waterman_matrix[i][j]) {
                sources[i][j] |= static_cast<int>(align_sources::match_align);
            }
        }
    }

    // todo: backtrace
}

void alignment::print_smith_waterman() {
    if (smith_waterman_matrix.empty()) {
        smith_waterman();
    }
    std::cout << "Smith-Waterman matrix:" << std::endl;
    print(smith_waterman_matrix);
}

void alignment::repeated_local_alignment(const int threshold) {
    if (!repeated_local_alignment_matrix.empty()) {
        repeated_local_alignment_matrix.clear();
    } else {
        repeated_local_alignment_matrix.resize(m + 1, std::vector(n + 1, 0));
    }
    std::vector sources(m + 1, std::vector(n + 1, static_cast<int>(align_sources::threshold_align)));
    int prev_max = 0;
    for (int j = 1; j <= n; ++j) {
        sources[0][j] = static_cast<int>(align_sources::threshold_align);
        repeated_local_alignment_matrix[0][j] = std::max(repeated_local_alignment_matrix[0][j-1],
                                                         prev_max - threshold);
        prev_max = 0;
        for (int i = 1; i <= m; ++i) {
            repeated_local_alignment_matrix[i][j] = repeated_local_alignment_matrix[0][j];

            if (const int _delete = repeated_local_alignment_matrix[i - 1][j] - d;
                _delete > repeated_local_alignment_matrix[i][j]) {
                repeated_local_alignment_matrix[i][j] = _delete;
                sources[i][j] = static_cast<int>(align_sources::delete_align);
            } else {
                sources[i][j] |= static_cast<int>(align_sources::delete_align);
            }

            if (const int _insert = repeated_local_alignment_matrix[i][j - 1] - d;
                _insert > repeated_local_alignment_matrix[i][j]) {
                repeated_local_alignment_matrix[i][j] = _insert;
                sources[i][j] = static_cast<int>(align_sources::insert_align);
            } else {
                sources[i][j] |= static_cast<int>(align_sources::insert_align);
            }

            if (const int _match = repeated_local_alignment_matrix[i - 1][j - 1] + similarity(
                                       seq1[j - 1], seq2[i - 1]);
                _match > repeated_local_alignment_matrix[i][j]) {
                repeated_local_alignment_matrix[i][j] = _match;
                sources[i][j] = static_cast<int>(align_sources::match_align);
            } else {
                sources[i][j] |= static_cast<int>(align_sources::match_align);
            }

            prev_max = std::max(prev_max, repeated_local_alignment_matrix[i][j]);
        }
    }

    // todo: backtrace
}

void alignment::print_repeated_local_alignment(const int threshold) {
    if (repeated_local_alignment_matrix.empty()) {
        repeated_local_alignment(threshold);
    }
    std::cout << "Repeated local alignment matrix:" << std::endl;
    print(repeated_local_alignment_matrix);
}

int alignment::wagner_fischer() {
    if (!wagner_fischer_matrix.empty()) {
        throw std::runtime_error("Wagner-Fischer matrix is already computed");
    }
    wagner_fischer_matrix.resize(m + 1, std::vector(n + 1, 0));
    for (int i = 1; i <= m; ++i) {
        wagner_fischer_matrix[i][0] = i;
    }
    for (int j = 1; j <= n; ++j) {
        wagner_fischer_matrix[0][j] = j;
    }
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            wagner_fischer_matrix[i][j] = std::min({
                wagner_fischer_matrix[i - 1][j] + 1,
                wagner_fischer_matrix[i][j - 1] + 1,
                wagner_fischer_matrix[i - 1][j - 1] +
                (seq1[j - 1] != seq2[i - 1])
            });
        }
    }
    return wagner_fischer_matrix[m][n];
}

void alignment::print_wagner_fischer() {
    if (wagner_fischer_matrix.empty()) {
        wagner_fischer();
    }
    std::cout << "Wagner-Fischer matrix:" << std::endl;
    print(wagner_fischer_matrix);
}

void alignment::print(const std::vector<std::vector<int> > &dp) {
    int width = 1;
    for (const auto &row: dp) {
        for (const auto &cell: row) {
            width = std::max(width, static_cast<int>(std::to_string(cell).size()));
        }
    }
    width++;
    std::cout << std::setw(width) << ' ' << std::setw(width) << ' ';
    for (const char c: seq1) {
        std::cout << std::setw(width) << c;
    }
    for (int i = 0; i <= m; ++i) {
        if (i == 0) {
            std::cout << std::endl << std::setw(width) << ' ';
        } else {
            std::cout << std::setw(width) << seq2[i - 1];
        }
        for (int j = 0; j <= n; ++j) {
            std::cout << std::setw(width) << dp[i][j];
        }
        std::cout << std::endl;
    }
}

int alignment::similarity(const char a, const char b) {
    int idx_a = 24, idx_b = 24;
    for (int i = 0; i < 24; i++) {
        if (Acids[i] == static_cast<char>(std::toupper(a))) {
            idx_a = i;
        }
        if (Acids[i] == static_cast<char>(std::toupper(b))) {
            idx_b = i;
        }
    }
    return ScoreMatrix[static_cast<int>(sm)][idx_a * 25 + idx_b];
}
