// differentials.h

#ifndef DIFFERENTIALS_H
#define DIFFERENTIALS_H

#include "../knotkit/knotkit.h"

template<class R>
class Differentials {
public:
    Differentials(const knot_diagram& kd, bool reduced = false);

    mod_map<R> compute_differential();

    void display_differential();

    void display_differential_component(int h, int q, int delta_q = 1);

private:
    knot_diagram kd_;
    bool reduced_;
    cube<R> c_;
    ptr<const module<R>> C_;
    mod_map<R> d_;
};

#endif // DIFFERENTIALS_H