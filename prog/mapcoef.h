// Copyright 2019 Global Phasing Ltd.
#pragma once

#include <gemmi/grid.hpp>      // for Grid
#include <optionparser.h>

// used by sf2map and blobs
enum MapOptions { Diff=4, Section, FLabel, PhLabel, WeightLabel, GridDims,
                  ExactDims, Sample, AxesZyx, GridQuery, TimingFft,
                  AfterMapOptions };

extern const option::Descriptor MapUsage[];

gemmi::Grid<float>
read_sf_and_fft_to_map(const char* input_path,
                       const std::vector<option::Option>& options,
                       FILE* output,
                       bool oversample_by_default=false);

// avoid including both mapcoef.h and monlib_opt.h
#ifdef GEMMI_OPTIONS_4
#error Conflicting headers
#endif
#define GEMMI_OPTIONS_4
