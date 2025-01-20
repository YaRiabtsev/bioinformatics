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

#include <benchmark/benchmark.h>
#include <random>
#include <string>
#include "findstr.hpp"

std::string generate_random_string(const size_t length) {
    const std::string characters = "ACGT";
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<> distribution(0, static_cast<int>(characters.size()) - 1);

    std::string random_string;
    for (size_t i = 0; i < length; ++i) {
        random_string += characters[distribution(generator)];
    }
    return random_string;
}

static void BM_KnuthMorrisPratt_Prefix(benchmark::State &state) {
    const size_t text_length = state.range(0);
    const size_t pattern_length = state.range(1);
    const std::string text = generate_random_string(text_length);
    const std::string pattern = generate_random_string(pattern_length);
    for (auto _: state) {
        std::vector<int> result;
        knuth_morris_pratt(text, pattern, result);
    }
}

static void BM_KnuthMorrisPratt_Z(benchmark::State &state) {
    const size_t text_length = state.range(0);
    const size_t pattern_length = state.range(1);
    const std::string text = generate_random_string(text_length);
    const std::string pattern = generate_random_string(pattern_length);
    for (auto _: state) {
        std::vector<int> result;
        knuth_morris_pratt(text, pattern, result, false);
    }
}

static void BM_BoyerMoore(benchmark::State &state) {
    const size_t text_length = state.range(0);
    const size_t pattern_length = state.range(1);
    const std::string text = generate_random_string(text_length);
    const std::string pattern = generate_random_string(pattern_length);
    for (auto _: state) {
        std::vector<int> result;
        boyer_moore(text, pattern, result);
    }
}

static void BM_Search(benchmark::State &state) {
    const size_t text_length = state.range(0);
    const size_t pattern_length = state.range(1);
    const std::string text = generate_random_string(text_length);
    const std::string pattern = generate_random_string(pattern_length);
    for (auto _: state) {
        std::vector<int> result;
        search(text, pattern, result);
    }
}

static void BM_StdFind(benchmark::State &state) {
    const size_t text_length = state.range(0);
    const size_t pattern_length = state.range(1);
    const std::string text = generate_random_string(text_length);
    const std::string pattern = generate_random_string(pattern_length);
    for (auto _: state) {
        std::vector<int> result;
        std_find(text, pattern, result);
    }
}


void CustomArguments(benchmark::internal::Benchmark *b) {
    for (int i = 10, j = 1000, k = 0; k < 3; ++k, i *= 10, j *= 1000) {
        b->Args({j, i});
    }
}

BENCHMARK(BM_KnuthMorrisPratt_Prefix)->Apply(CustomArguments);
BENCHMARK(BM_KnuthMorrisPratt_Z)->Apply(CustomArguments);
BENCHMARK(BM_BoyerMoore)->Apply(CustomArguments);
BENCHMARK(BM_Search)->Apply(CustomArguments);
BENCHMARK(BM_StdFind)->Apply(CustomArguments);

BENCHMARK_MAIN();
