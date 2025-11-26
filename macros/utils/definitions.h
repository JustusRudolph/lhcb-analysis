#pragma once

#include <string>


namespace Utils::Definitions {
  std::string stackRoot = std::getenv("STACK_ROOT");
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

  constexpr float pi_f_float = 3.141592654f;
  constexpr float c_f_float = 299792458.f;  // speed of light in m/s
  constexpr float c_f_float_mmns = c_f_float * 1e3f * 1e-9f;  // speed of light in mm/ns
  constexpr float inv_c_mmns = 1.f / c_f_float_mmns;

  constexpr unsigned kModules = 64;
}