// Just ported functions from Allen to accurately "reconstruct" the same values as would be
// seen in the running version of Allen
#pragma once

#include <cmath>

#include "definitions.h"
#include "Hit.h"

namespace Utils {
  namespace Functions {

    int hit_phi_float_to_16(const float hit_phi) {
      float max_output_value_i16 = 1 << 15;
      float convert_factor_i16 = max_output_value_i16 / Definitions::pi_f_float;
      float hit_phi_i16_float = hit_phi * convert_factor_i16;
      
      return static_cast<int16_t>(hit_phi_i16_float);;
    }

    // return estimated z of the next module
    float get_next_module_z(const Hit::BaseHit h0,
                             const Hit::BaseHit h1,
                             const std::vector<double>& moduleToZ) {
      // check which direction track is in
      unsigned h1_module_number = (h1.id >> 12) & 0x3F;
      if (h1.z < h0.z && h1_module_number > 1) {  // moving backwards (as we should in normal case)
        return moduleToZ[h1_module_number - 2];
      } else if (h1.z > h0.z && h1_module_number < 62) {  // moving forwards (can't be at edge (62&63))
        return moduleToZ[h1_module_number + 2];
      } else {
        return 1e6f;  // no next module, out of scope (800mm max rn)
      }
    }

    // return estimated absolute time difference to next module/layer
    float estimate_abs_dt(const Hit::BaseHit h0,
                          const Hit::BaseHit h1,
                          const float next_module_z,
                          unsigned order=0) {
      if (next_module_z > 1000.) {
        return 0.;  // no next module, so expected dt is zero
      }

      // at first order we use x and y too
      float z_mid = (h1.z + h0.z) / 2;
      float dz_next_module = next_module_z - z_mid;
      float dt;
      if (order == 0) {
        dt = std::abs(dz_next_module / Definitions::c_f_float_mmns);
      } else if (order == 1) {
        float dx = h1.x - h0.x;
        float dy = h1.y - h0.y;
        float dz = h1.z - h0.z;
        float dx_next_module = (dx/dz) * dz_next_module;
        float dy_next_module = (dy/dz) * dz_next_module;
        float dr_sq = dx_next_module * dx_next_module +
                      dy_next_module * dy_next_module +
                      dz_next_module * dz_next_module;
        float dt_sq = dr_sq / (Definitions::c_f_float_mmns * 
                               Definitions::c_f_float_mmns);
        dt = std::sqrt(dt_sq);
      } else { // add more orders later (like how far back to go in propagation)
        throw std::invalid_argument("Order must be 0 or 1");
      }
      return dt;
    }

    // return estimated x, y and thus phi of the next hit with module to Z vector
    std::tuple<float, float, float> estimateNextPhi(const Hit::BaseHit h0,
                                                    const Hit::BaseHit h1,
                                                    const float next_module_z) {
      if (next_module_z > 1000.) {
        return {0., 0., 0.};  // no next module, put expected position in beamline
      }
      float dz_to_next_module = next_module_z - h0.z;
      
      float dz = h1.z - h0.z;
      float dx = (h1.x - h0.x);
      float dy = (h1.y - h0.y);
      float dxdz = dx / dz;
      float dydz = dy / dz;
      // take difference to h0, that way we get longer interpolation (same as in real alg)
      float predx = dxdz * dz_to_next_module;
      float predy = dydz * dz_to_next_module;
      float x_prediction = h0.x + predx;
      float y_prediction = h0.y + predy;
      float track_extrapolation_phi = TMath::ATan2(y_prediction, x_prediction);

      return {x_prediction, y_prediction, track_extrapolation_phi};
    }
  }
}