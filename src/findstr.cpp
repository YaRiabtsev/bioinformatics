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

#include <unordered_map>
#include "findstr.hpp"


void knuth_morris_pratt(const std::string &text, const std::string &pattern, std::vector<int> &answer,
                        const bool prefix) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pattern.size());
    if (m == 0) {
        answer.push_back(-1);
        return;
    }
    std::vector<int> dp;
    if (prefix) {
        prefix_function(pattern + "#" + text, dp);
    } else {
        z_function(pattern + "#" + text, dp);
    }
    for (int i = 0; i < n; i++) {
        if (dp[m + i + 1] == m) {
            if (prefix) {
                answer.push_back(i - m + 1);
            } else {
                answer.push_back(i);
            }
        }
    }
    if (answer.empty()) {
        answer.push_back(-1);
    }
}

bool is_prefix(const std::string &pattern, const int p) {
    int j = 0;
    for (int i = p; i < pattern.size(); i++) {
        if (pattern[i] != pattern[j]) {
            return false;
        }
        j++;
    }
    return true;
}

int suffix_length(const std::string &pattern, const int p) {
    int length = 0;
    int i = p;
    int j = static_cast<int>(pattern.size()) - 1;
    while (i >= 0 && pattern[i] == pattern[j]) {
        length++;
        i--;
        j--;
    }
    return length;
}

void boyer_moore(const std::string &text, const std::string &pattern, std::vector<int> &answer) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pattern.size());
    if (m == 0) {
        answer.push_back(-1);
        return;
    }
    std::unordered_map<char, int> bad_char;
    bad_char[pattern[m - 1]] = m;
    for (int i = 0; i < m - 1; i++) {
        bad_char[pattern[i]] = m - 1 - i;
    }
    std::vector<int> good_suffix(m);
    int last_prefix = m;
    for (int i = m - 1; i >= 0; i--) {
        if (is_prefix(pattern, i + 1)) {
            last_prefix = i + 1;
        }
        good_suffix[m - 1 - i] = last_prefix - i + m - 1;
    }
    for (int i = 0; i < m - 1; i++) {
        const int suffix_len = suffix_length(pattern, i);
        good_suffix[suffix_len] = m - 1 - i + suffix_len;
    }
    for (int i = m - 1; i < n;) {
        int j = m - 1;
        while (j >= 0 && text[i] == pattern[j]) {
            j--;
            i--;
        }
        if (j < 0) {
            answer.push_back(i + 1);
            i += m + 1;
        } else {
            i += std::max(good_suffix[m - 1 - j], bad_char[text[j]]);
        }
    }
    if (answer.empty()) {
        answer.push_back(-1);
    }
}

void search(const std::string &text, const std::string &pattern, std::vector<int> &answer) {
    const int n = static_cast<int>(text.size());
    const int m = static_cast<int>(pattern.size());
    if (m == 0) {
        answer.push_back(-1);
        return;
    }
    for (int offset = 0; offset < n - m + 1; offset++) {
        bool found = true;
        for (int i = 0; i < m; i++) {
            if (text[offset + i] != pattern[i]) {
                found = false;
                break;
            }
        }
        if (found) {
            answer.push_back(offset);
        }
    }
    if (answer.empty()) {
        answer.push_back(-1);
    }
}

void std_find(const std::string &text, const std::string &pattern, std::vector<int> &answer) {
    if (const int m = static_cast<int>(pattern.size()); m == 0) {
        answer.push_back(-1);
        return;
    }
    int idx = -1;
    while (true) {
        idx = static_cast<int>(text.find(pattern, idx + 1));
        if (idx == std::string::npos) {
            break;
        }
        answer.push_back(idx);
    }
    if (answer.empty()) {
        answer.push_back(-1);
    }
}

void prefix_function(const std::string &test, std::vector<int> &dp) {
    const int n = static_cast<int>(test.size());
    dp.resize(n);
    dp[0] = 0;
    for (int i = 1; i < n; i++) {
        int k = dp[i - 1];
        while (k > 0 && test[i] != test[k]) {
            k = dp[k - 1];
        }
        if (test[i] == test[k]) {
            k++;
        }
        dp[i] = k;
    }
}

void z_function(const std::string &test, std::vector<int> &dp) {
    const int n = static_cast<int>(test.size());
    dp.resize(n);
    dp[0] = 0;
    int left = 0, right = 0;
    for (int i = 1; i < n; i++) {
        dp[i] = std::max(0, std::min(dp[i - left], right - i));
        while (i + dp[i] < n && test[dp[i]] == test[i + dp[i]]) {
            dp[i]++;
        }
        if (i + dp[i] > right) {
            left = i;
            right = i + dp[i];
        }
    }
}
