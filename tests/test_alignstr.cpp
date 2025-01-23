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

#include <gtest/gtest.h>
#include "alignstr.hpp"

TEST(AlignmentStringTest, BaseTest) {
    alignment align("kitten", "sitting");
    EXPECT_EQ(align.wagner_fischer(), 3);

    align = alignment("Saturday", "Sunday");
    EXPECT_EQ(align.wagner_fischer(), 3);
}

TEST(AlignmentStringPrintTest, PrintTest) {
    alignment align("HEAGAWGHEE", "PAWHEAE", score_matrix::blosum50, 8);

    std::vector<std::pair<std::string, std::string> > best_align_pairs;
    align.needleman_wunsch(best_align_pairs);
    const std::vector<std::pair<std::string, std::string> > best_align_pairs_expected = {
        {"HEAGAWGHE-E", "--P-AW-HEAE"},
        {"HEAGAWGHE-E", "-P--AW-HEAE"},
        {"HEAGAWGHE-E", "-PA--W-HEAE"}
    };
    EXPECT_EQ(best_align_pairs, best_align_pairs_expected);

    std::vector<std::vector<std::pair<std::string, std::string> > > best_local_align_pairs;
    EXPECT_EQ(align.smith_waterman(best_local_align_pairs), 28);
    const std::vector<std::vector<std::pair<std::string, std::string> > > best_local_align_pairs_expected = {
        {
            {"AWGHE", "AW-HE"}
        }
    };
    EXPECT_EQ(best_local_align_pairs, best_local_align_pairs_expected);

    std::vector<std::pair<std::string, std::string> > best_rlocal_align_pairs;
    align.repeated_local_alignment(best_rlocal_align_pairs);
    const std::vector<std::pair<std::string, std::string> > best_rlocal_align_pairs_expected = {
        {"HEAGAWGHEE", "HEA_AW-HE_"}
    };
    EXPECT_EQ(best_rlocal_align_pairs, best_rlocal_align_pairs_expected);
}
