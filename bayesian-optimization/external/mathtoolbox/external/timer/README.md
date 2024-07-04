# timer

![GitHub](https://img.shields.io/github/license/yuki-koyama/timer)

Execution timer for measuring elapsed times for code blocks written in C++11

## Usage

Suppose that you want to measure the execution time for a function in a C++ code:
```
some_function();
```
By using `timer`, this is easily achieved:
```
{
  timer::Timer t;
  some_function();
}
```
It will automatically display the elapsed time in `std::cout`, for example,
```
[ timer : 24.601 ms ]
```

## Installation

`timer` is a header-only, single-file library. It can be used by just copying `timer.hpp` and pasting it into your project.

Alternatively, it can be installed using `cmake`. If your project is also managed using `cmake`, `ExternalProject` commands are useful for including `timer` to your project. Another option is to use the `add_subdirectory` command to include this repository.

## Licensing

MIT License.
