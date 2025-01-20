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

#ifndef FINDSTR_HPP
#define FINDSTR_HPP

#include <string>
#include <vector>

// todo: Aho–Corasick algorithm
// todo: Rabin-Karp algorithm
// todo: Colussi algorithm

void knuth_morris_pratt(const std::string &text, const std::string &pattern, std::vector<int>& result,
                        bool prefix = true);

void boyer_moore(const std::string &text, const std::string &pattern, std::vector<int> &result);

void search(const std::string &text, const std::string &pattern, std::vector<int> &result);

void std_find(const std::string &text, const std::string &pattern, std::vector<int> &result);

void prefix_function(const std::string &test, std::vector<int> &dp);

void z_function(const std::string &test, std::vector<int> &dp);

#endif //FINDSTR_HPP
