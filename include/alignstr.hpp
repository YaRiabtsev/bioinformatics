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

#ifndef ALIGNSTR_HPP
#define ALIGNSTR_HPP
#include "tables.hpp"

class alignment {
public:
    explicit alignment(std::string seq_a, std::string seq_b, score_matrix sm = score_matrix::blosum50, int d = 6);


    void needleman_wunsch(std::vector<std::pair<std::string, std::string> > &best_align_pairs);

    void print_needleman_wunsch();

    int smith_waterman(std::vector<std::vector<std::pair<std::string, std::string> > > &best_align_pairs);

    void print_smith_waterman();

    void repeated_local_alignment(std::vector<std::pair<std::string, std::string> > &best_align_pairs,
                                  int threshold = 20);

    void print_repeated_local_alignment();

    int wagner_fischer();

    void print_wagner_fischer();

private:
    int d; /// gap penalty score
    score_matrix sm;
    std::string seq1, seq2;
    size_t n, m;
    std::vector<std::vector<int> > needleman_wunsch_matrix;
    std::vector<std::vector<int> > smith_waterman_matrix;
    std::vector<std::vector<int> > repeated_local_alignment_matrix;
    std::vector<std::vector<int> > wagner_fischer_matrix;


    enum class align_sources : int {
        null_align = 0, match_align = 1, delete_align = 2, insert_align = 4, threshold_align = 8
    };

    void print(const std::vector<std::vector<int> > &dp);

    int similarity(char a, char b);

    void backtrace(const std::vector<std::vector<int> > &sources, size_t i, size_t j,
                   std::vector<std::pair<std::string, std::string> > &pairs);
};

#endif //ALIGNSTR_HPP
