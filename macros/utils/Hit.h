#pragma once

#include <TMath.h>

namespace Hit {
  struct BaseHit {
    unsigned id;
    float x, y, z, t;
    // default constructor
    BaseHit() : id(0), x(0), y(0), z(0), t(0) {}
    // constructor with parameters
    BaseHit(unsigned id, float x, float y, float z, float t)
        : id(id), x(x), y(y), z(z), t(t) {}
    // copy constructor
    BaseHit(const BaseHit& other) : id(other.id), x(other.x), y(other.y), z(other.z), t(other.t) {}
  };

  float estimateNextPhi(BaseHit h0, BaseHit h1, float nextZ) {
    float dz = h1.z - h0.z;
    float dx = (h1.x - h0.x);
    float dy = (h1.y - h0.y);
    float dxdz = dx / dz;
    float dydz = dy / dz;
    // take difference to h0, that way we get longer interpolation (same as in real alg)
    float dz_to_next_module = nextZ - h0.z;
    float predx = dxdz * dz_to_next_module;
    float predy = dydz * dz_to_next_module;
    float x_prediction = h0.x + predx;
    float y_prediction = h0.y + predy;
    float track_extrapolation_phi = TMath::ATan2(y_prediction, x_prediction);

    return track_extrapolation_phi;
  }
}