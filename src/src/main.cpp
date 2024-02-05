// MIT License

// Copyright (c) 2021 Daniel Kowalczyk

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <docopt/docopt.h>  // for DocoptArgumentError, DocoptLanguageError
#include <fmt/core.h>       // for print
#include <cstdio>           // for stderr
#include <exception>        // for exception
#include "Problem.h"        // for Problem
#include "Usage.hpp"        // for USAGE

auto main(int argc, const char** argv) -> int {
    int val = 0;
    try {
        Problem problem(argc, argv);
    } catch (docopt::DocoptExitHelp const&) {
        fmt::print("{}", USAGE);
    } catch (docopt::DocoptExitVersion const&) {
        fmt::print("PM 0.1\n");
    } catch (docopt::DocoptLanguageError const& error) {
        fmt::print(stderr, "Docopt usage string could not be parsed\n");
        fmt::print(stderr, "{}\n", error.what());
    } catch (docopt::DocoptArgumentError const& error) {
        fmt::print("{}\n", error.what());
        fmt::print("{}", USAGE);
    } catch (std::exception& e) {
        fmt::print(stderr, "{}\n", e.what());
    } catch (...) {
        fmt::print(stderr, "error: unknown exceptions\n");
    }
    return val;
}
