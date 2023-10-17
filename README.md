# lineIntSquaredExponential

Comparison code for the paper at https://arxiv.org/abs/1812.07319

## Dependencies

Install [gsl](https://www.gnu.org/software/gsl/doc/html/intro.html).

If building on MAC OS

```bash
brew install gsl
```

If building on Ubuntu

```bash
sudo apt-get install gsl
```

## Installation

```bash
cmake -G Ninja -B build
cmake --build build
```

## Usage

Two example Matlab scripts are included

1. `comparison_set_2.m`: Runs set two from the paper and plots the histogram
2. `comparison_all_sets.m`: Runs all sets from the paper and plots the output. Note, the log10 of the error is being plotted which for set 4 is -inf and so displays as a blank point on the plot
