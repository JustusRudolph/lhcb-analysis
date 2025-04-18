#pragma once

#include <string>


namespace Utils::Definitions {
  std::string stackRoot = std::getenv("STACK_ROOT");
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");
}