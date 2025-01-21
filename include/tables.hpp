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

#ifndef TABLES_HPP
#define TABLES_HPP

#include <string>

enum class score_matrix : int {
    blosum45 = 0, blosum50, blosum62, blosum80, blosum90, pam250
};

constexpr char Acids[26] = "ARNDCQEGHILKMFPSTWYVBJZX*";

constexpr int ScoreMatrix[6][25 * 25] = {
    {
        /// blosum45:
        5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2, -2, 0, -1, -1, -1, -1, -5, -2, 7, 0, -1, -3,
        1, 0, -2, 0, -3, -2, 3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3, 1, -1, -5, -1, 0, 6, 2, -2, 0, 0, 0, 1, -2, -3,
        0, -2, -2, -2, 1, 0, -4, -2, -3, 5, -3, 0, -1, -5, -2, -1, 2, 7, -3, 0, 2, -1, 0, -4, -3, 0, -3, -4, -1, 0, -1,
        -4, -2, -3, 6, -3, 1, -1, -5, -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1,
        -2, -2, -3, -1, -5, -1, 1, 0, 0, -3, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3, 0, -2, 4, -1, -5, -1,
        0, 0, 2, -3, 2, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3, 1, -3, 5, -1, -5, 0, -2, 0, -1, -3, -2, -2,
        7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -5, -2, 0, 1, 0, -3, 1, 0, -2, 10, -3, -2, -1,
        0, -2, -2, -1, -2, -3, 2, -3, 0, -2, 0, -1, -5, -1, -3, -2, -4, -3, -2, -3, -4, -3, 5, 2, -3, 2, 0, -2, -2, -1,
        -2, 0, 3, -3, 4, -3, -1, -5, -1, -2, -3, -3, -2, -2, -2, -3, -2, 2, 5, -3, 2, 1, -3, -3, -1, -2, 0, 1, -3, 4,
        -2, -1, -5, -1, 3, 0, 0, -3, 1, 1, -2, -1, -3, -3, 5, -1, -3, -1, -1, -1, -2, -1, -2, 0, -3, 1, -1, -5, -1, -1,
        -2, -3, -2, 0, -2, -2, 0, 2, 2, -1, 6, 0, -2, -2, -1, -2, 0, 1, -2, 2, -1, -1, -5, -2, -2, -2, -4, -2, -4, -3,
        -3, -2, 0, 1, -3, 0, 8, -3, -2, -1, 1, 3, 0, -3, 1, -3, -1, -5, -1, -2, -2, -1, -4, -1, 0, -2, -2, -2, -3, -1,
        -2, -3, 9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5, 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -3, -1, -2, -2, -1, 4, 2,
        -4, -2, -1, 0, -2, 0, -1, -5, 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1, 2, 5, -3, -1, 0, 0, -1,
        -1, -1, -5, -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2, 1, -3, -4, -3, 15, 3, -3, -4, -2, -2, -1, -5,
        -2, -1, -2, -2, -3, -1, -2, -3, 2, 0, 0, -1, 0, 3, -3, -2, -1, 3, 8, -1, -2, 0, -2, -1, -5, 0, -2, -3, -3, -1,
        -3, -3, -3, -3, 3, 1, -2, 1, 0, -3, -1, 0, -3, -1, 5, -3, 2, -3, -1, -5, -1, -1, 5, 6, -2, 0, 1, -1, 0, -3, -3,
        0, -2, -3, -2, 0, 0, -4, -2, -3, 5, -3, 1, -1, -5, -1, -3, -3, -3, -2, -2, -3, -4, -2, 4, 4, -3, 2, 1, -3, -2,
        -1, -2, 0, 2, -3, 4, -2, -1, -5, -1, 1, 0, 1, -3, 4, 5, -2, 0, -3, -2, 1, -1, -3, -1, 0, -1, -2, -2, -3, 1, -2,
        5, -1, -5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5,
        -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1
    },
    {
        /// blosum50:
        5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -5, -2, 7, -1, -2,
        -4, 1, 0, -3, 0, -4, -3, 3, -2, -3, -3, -1, -1, -3, -1, -3, -1, -3, 0, -1, -5, -1, -1, 7, 2, -2, 0, 0, 0, 1, -3,
        -4, 0, -2, -4, -2, 1, 0, -4, -2, -3, 5, -4, 0, -1, -5, -2, -2, 2, 8, -4, 0, 2, -1, -1, -4, -4, -1, -4, -5, -1,
        0, -1, -5, -3, -4, 6, -4, 1, -1, -5, -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3,
        -1, -3, -2, -3, -1, -5, -1, 1, 0, 0, -3, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3, 0, -3, 4, -1, -5,
        -1, 0, 0, 2, -3, 2, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3, 1, -3, 5, -1, -5, 0, -3, 0, -1, -3, -2,
        -3, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4, -1, -4, -2, -1, -5, -2, 0, 1, -1, -3, 1, 0, -2, 10, -4,
        -3, 0, -1, -1, -2, -1, -2, -3, 2, -4, 0, -3, 0, -1, -5, -1, -4, -3, -4, -2, -3, -4, -4, -4, 5, 2, -3, 2, 0, -3,
        -3, -1, -3, -1, 4, -4, 4, -3, -1, -5, -2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5, -3, 3, 1, -4, -3, -1, -2, -1, 1,
        -4, 4, -3, -1, -5, -1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6, -2, -4, -1, 0, -1, -3, -2, -3, 0, -3, 1, -1, -5,
        -1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7, 0, -3, -2, -1, -1, 0, 1, -3, 2, -1, -1, -5, -3, -3, -4, -5, -2,
        -4, -3, -4, -1, 0, 1, -4, 0, 8, -4, -3, -2, 1, 4, -1, -4, 1, -4, -1, -5, -1, -3, -2, -1, -4, -1, -1, -2, -2, -3,
        -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5, 1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3,
        -1, 5, 2, -4, -2, -2, 0, -3, 0, -1, -5, 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5, -3, -2,
        0, 0, -1, -1, -1, -5, -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15, 2, -3, -5, -2, -2,
        -1, -5, -2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8, -1, -3, -1, -2, -1, -5, 0, -3,
        -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5, -3, 2, -3, -1, -5, -2, -1, 5, 6, -3, 0, 1,
        -1, 0, -4, -4, 0, -3, -4, -2, 0, 0, -5, -3, -3, 6, -4, 1, -1, -5, -2, -3, -4, -4, -2, -3, -3, -4, -3, 4, 4, -3,
        2, 1, -3, -3, -1, -2, -1, 2, -4, 4, -3, -1, -5, -1, 0, 0, 1, -3, 4, 5, -2, 0, -3, -3, 1, -1, -4, -1, 0, -1, -2,
        -2, -3, 1, -3, 5, -1, -5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
        1
    },
    {
        /// blosum62:
        4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1, -1, -4, -1, 5, 0, -2, -3,
        1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2, 0, -1, -4, -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3,
        0, -2, -3, -2, 1, 0, -4, -2, -3, 4, -3, 0, -1, -4, -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0,
        -1, -4, -3, -3, 4, -3, 1, -1, -4, 0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1,
        -3, -1, -3, -1, -4, -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, -2, 4, -1, -4, -1,
        0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, -3, 4, -1, -4, 0, -2, 0, -1, -3, -2, -2,
        6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -4, -2, -1, -4, -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1,
        -2, -1, -2, -1, -2, -2, 2, -3, 0, -3, 0, -1, -4, -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1,
        -3, -1, 3, -3, 3, -3, -1, -4, -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, 3,
        -3, -1, -4, -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, -3, 1, -1, -4, -1, -1,
        -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, 2, -1, -1, -4, -2, -3, -3, -3, -2, -3, -3,
        -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, 0, -3, -1, -4, -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1,
        -2, -4, 7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4, 1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1,
        -3, -2, -2, 0, -2, 0, -1, -4, 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1,
        -1, -1, -4, -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -2, -2, -1, -4,
        -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -1, -2, -1, -4, 0, -3, -3, -3,
        -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, 2, -2, -1, -4, -2, -1, 4, 4, -3, 0, 1, -1, 0, -3,
        -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, -3, 0, -1, -4, -1, -2, -3, -3, -1, -2, -3, -4, -3, 3, 3, -3, 2, 0, -3,
        -2, -1, -2, -1, 2, -3, 3, -3, -1, -4, -1, 0, 0, 1, -3, 4, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -2, -2, -2, 0,
        -3, 4, -1, -4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1
    },
    {
        /// blosum80:
        5, -2, -2, -2, -1, -1, -1, 0, -2, -2, -2, -1, -1, -3, -1, 1, 0, -3, -2, 0, -2, -2, -1, -1, -6, -2, 6, -1, -2,
        -4, 1, -1, -3, 0, -3, -3, 2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3, 0, -1, -6, -2, -1, 6, 1, -3, 0, -1, -1, 0,
        -4, -4, 0, -3, -4, -3, 0, 0, -4, -3, -4, 5, -4, 0, -1, -6, -2, -2, 1, 6, -4, -1, 1, -2, -2, -4, -5, -1, -4, -4,
        -2, -1, -1, -6, -4, -4, 5, -5, 1, -1, -6, -1, -4, -3, -4, 9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3,
        -3, -1, -4, -2, -4, -1, -6, -1, 1, 0, -1, -4, 6, 2, -2, 1, -3, -3, 1, 0, -4, -2, 0, -1, -3, -2, -3, 0, -3, 4,
        -1, -6, -1, -1, -1, 1, -5, 2, 6, -3, 0, -4, -4, 1, -2, -4, -2, 0, -1, -4, -3, -3, 1, -4, 5, -1, -6, 0, -3, -1,
        -2, -4, -2, -3, 6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6, -2, 0, 0, -2, -4, 1, 0,
        -3, 8, -4, -3, -1, -2, -2, -3, -1, -2, -3, 2, -4, -1, -4, 0, -1, -6, -2, -3, -4, -4, -2, -3, -4, -5, -4, 5, 1,
        -3, 1, -1, -4, -3, -1, -3, -2, 3, -4, 3, -4, -1, -6, -2, -3, -4, -5, -2, -3, -4, -4, -3, 1, 4, -3, 2, 0, -3, -3,
        -2, -2, -2, 1, -4, 3, -3, -1, -6, -1, 2, 0, -1, -4, 1, 1, -2, -1, -3, -3, 5, -2, -4, -1, -1, -1, -4, -3, -3, -1,
        -3, 1, -1, -6, -1, -2, -3, -4, -2, 0, -2, -4, -2, 1, 2, -2, 6, 0, -3, -2, -1, -2, -2, 1, -3, 2, -1, -1, -6, -3,
        -4, -4, -4, -3, -4, -4, -4, -2, -1, 0, -4, 0, 6, -4, -3, -2, 0, 3, -1, -4, 0, -4, -1, -6, -1, -2, -3, -2, -4,
        -2, -2, -3, -3, -4, -3, -1, -3, -4, 8, -1, -2, -5, -4, -3, -2, -4, -2, -1, -6, 1, -1, 0, -1, -2, 0, 0, -1, -1,
        -3, -3, -1, -2, -3, -1, 5, 1, -4, -2, -2, 0, -3, 0, -1, -6, 0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1,
        -2, -2, 1, 5, -4, -2, 0, -1, -1, -1, -1, -6, -3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2, 0, -5, -4, -4,
        11, 2, -3, -5, -3, -3, -1, -6, -2, -3, -3, -4, -3, -2, -3, -4, 2, -2, -2, -3, -2, 3, -4, -2, -2, 2, 7, -2, -3,
        -2, -3, -1, -6, 0, -3, -4, -4, -1, -3, -3, -4, -4, 3, 1, -3, 1, -1, -3, -2, 0, -3, -2, 4, -4, 2, -3, -1, -6, -2,
        -1, 5, 5, -4, 0, 1, -1, -1, -4, -4, -1, -3, -4, -2, 0, -1, -5, -3, -4, 5, -4, 0, -1, -6, -2, -3, -4, -5, -2, -3,
        -4, -5, -4, 3, 3, -3, 2, 0, -4, -3, -1, -3, -2, 2, -4, 3, -3, -1, -6, -1, 0, 0, 1, -4, 4, 5, -3, 0, -4, -3, 1,
        -1, -4, -2, 0, -1, -3, -3, -3, 0, -3, 5, -1, -6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
        -6, -6, -6, -6, -6, 1
    },
    {
        /// blosum90:
        5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4, -3, -1, -2, -2, -1, -1, -6, -2, 6, -1, -3,
        -5, 1, -1, -3, 0, -4, -3, 2, -2, -4, -3, -1, -2, -4, -3, -3, -2, -3, 0, -1, -6, -2, -1, 7, 1, -4, 0, -1, -1, 0,
        -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, 5, -4, -1, -1, -6, -3, -3, 1, 7, -5, -1, 1, -2, -2, -5, -5, -1, -4, -5,
        -3, -1, -2, -6, -4, -5, 5, -5, 1, -1, -6, -1, -5, -4, -5, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4,
        -4, -2, -4, -2, -5, -1, -6, -1, 1, 0, -1, -4, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3, -3, -1, -3, 5,
        -1, -6, -1, -1, -1, 1, -6, 2, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5, -4, -3, 1, -4, 5, -1, -6, 0, -3, -1,
        -2, -4, -3, -3, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -5, -3, -1, -6, -2, 0, 0, -2, -5, 1, -1,
        -3, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1, -4, -1, -4, 0, -1, -6, -2, -4, -4, -5, -2, -4, -4, -5, -4, 5, 1,
        -4, 1, -1, -4, -3, -1, -4, -2, 3, -5, 3, -4, -1, -6, -2, -3, -4, -5, -2, -3, -4, -5, -4, 1, 5, -3, 2, 0, -4, -3,
        -2, -3, -2, 0, -5, 4, -4, -1, -6, -1, 2, 0, -1, -4, 1, 0, -2, -1, -4, -3, 6, -2, -4, -2, -1, -1, -5, -3, -3, -1,
        -3, 1, -1, -6, -2, -2, -3, -4, -2, 0, -3, -4, -3, 1, 2, -2, 7, -1, -3, -2, -1, -2, -2, 0, -4, 2, -2, -1, -6, -3,
        -4, -4, -5, -3, -4, -5, -5, -2, -1, 0, -4, -1, 7, -4, -3, -3, 0, 3, -2, -4, 0, -4, -1, -6, -1, -3, -3, -3, -4,
        -2, -2, -3, -3, -4, -4, -2, -3, -4, 8, -2, -2, -5, -4, -3, -3, -4, -2, -1, -6, 1, -1, 0, -1, -2, -1, -1, -1, -2,
        -3, -3, -1, -2, -3, -2, 5, 1, -4, -3, -2, 0, -3, -1, -1, -6, 0, -2, 0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1,
        -3, -2, 1, 6, -4, -2, -1, -1, -2, -1, -1, -6, -4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2, 0, -5, -4, -4,
        11, 2, -3, -6, -3, -4, -1, -6, -3, -3, -3, -4, -4, -3, -4, -5, 1, -2, -2, -3, -2, 3, -4, -3, -2, 2, 8, -3, -4,
        -2, -3, -1, -6, -1, -3, -4, -5, -2, -3, -3, -5, -4, 3, 0, -3, 0, -2, -3, -2, -1, -3, -3, 5, -4, 1, -3, -1, -6,
        -2, -2, 5, 5, -4, -1, 1, -2, -1, -5, -5, -1, -4, -4, -3, 0, -1, -6, -4, -4, 5, -5, 0, -1, -6, -2, -3, -4, -5,
        -2, -3, -4, -5, -4, 3, 4, -3, 2, 0, -4, -3, -2, -3, -2, 1, -5, 4, -4, -1, -6, -1, 0, -1, 1, -5, 5, 5, -3, 0, -4,
        -4, 1, -2, -4, -2, -1, -1, -4, -3, -3, 0, -4, 5, -1, -6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,
        -6, -6, -6, -6, -6, -6, -6, 1
    },
    {
        /// pam250:
        2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3, 0, 0, -1, 0, -1, -8, -2, 6, 0, -1, -4, 1, -1,
        -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, -1, -3, 0, -1, -8, 0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -3,
        0, 1, 0, -4, -2, -2, 2, -3, 1, -1, -8, 0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, 3,
        -3, 3, -1, -8, -2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -5, -1, -8,
        0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, -2, 3, -1, -8, 0, -1, 1, 3, -5, 2, 4,
        0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, -3, 3, -1, -8, 1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3,
        -5, 0, 1, 0, -7, -5, -1, 0, -4, 0, -1, -8, -1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0,
        -2, 1, -2, 2, -1, -8, -1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, -2, 3, -2, -1,
        -8, -2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, -3, 5, -3, -1, -8, -1, 3, 1, 0,
        -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, -3, 0, -1, -8, -1, 0, -2, -3, -5, -1, -2, -3, -2, 2,
        4, 0, 6, 0, -2, -2, -1, -4, -2, 2, -2, 3, -2, -1, -8, -3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5,
        -3, -3, 0, 7, -1, -4, 2, -5, -1, -8, 1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, -1,
        -2, 0, -1, -8, 1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2, 1, -2, -3, -1, 0, -2, 0, -1, -8, 1, -1, 0,
        0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3, -5, -3, 0, 0, -1, -1, -1, -8, -6, 2, -4, -7, -8, -5, -7, -7, -3,
        -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, -5, -3, -6, -1, -8, -3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2,
        7, -5, -3, -3, 0, 10, -2, -3, -1, -4, -1, -8, 0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6,
        -2, 4, -2, 2, -2, -1, -8, 0, -1, 2, 3, -4, 1, 3, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3, -3, 2, -1,
        -8, -1, -3, -3, -3, -5, -2, -3, -4, -2, 3, 5, -3, 3, 2, -2, -2, -1, -3, -1, 2, -3, 5, -2, -1, -8, 0, 0, 1, 3,
        -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, -2, 3, -1, -8, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,
        -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, 1
    }
};

#endif //TABLES_HPP
