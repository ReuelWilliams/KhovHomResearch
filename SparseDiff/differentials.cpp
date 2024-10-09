#include "../knotkit/knotkit.h"

template<class R>
class Differentials {
public:
    Differentials(const knot_diagram& kd, bool reduced = false);

    mod_map<R> compute_differential();

    void display_differential();

    void display_differential_component(int h, int q, int delta_q = 0);

private:
    knot_diagram kd_;
    bool reduced_;
    cube<R> c_;
    ptr<const module<R>> C_;
    mod_map<R> d_;
};

template<class R>
Differentials<R>::Differentials(const knot_diagram& kd, bool reduced)
    : kd_(kd), reduced_(reduced), c_(kd_, reduced_) {
    // Initialize the cube of resolutions
    C_ = c_.khC;
    d_ = c_.compute_d(1, 0, 0, 0, 0);
}

template<class R>
mod_map<R> Differentials<R>::compute_differential() {
    return d_;
}

template<class R>
void Differentials<R>::display_differential() {
    d_.display_self();
}

template<class R>
void Differentials<R>::display_differential_component(int h, int q, int delta_q) {
    // Extract the graded pieces
    grading from_grading(h, q);
    grading to_grading(h - 1, q + delta_q);

    // Get the submodules at the specified gradings
    auto C_from = C_->graded_piece(from_grading);
    auto C_to = C_->graded_piece(to_grading);

    if (C_from->dim() == 0 || C_to->dim() == 0) {
        std::cout << "No generators at the specified gradings.\n";
        return;
    }

    // Restrict the differential to these submodules
    mod_map<R> d_restricted = d_.restrict(C_from, C_to);

    // Display the restricted differential
    std::cout << "Differential from grading (" << h << ", " << q << ") "
              << "to (" << h - 1 << ", " << q + delta_q << "):\n";

    // Represent the differential as a matrix
    unsigned m = C_from->dim();
    unsigned n = C_to->dim();

    // Build the matrix
    std::vector<std::vector<R>> matrix(n, std::vector<R>(m, R(0)));

    for (unsigned i = 1; i <= m; ++i) {
        linear_combination<R> col = d_restricted.column(i);
        for (linear_combination_const_iter<R> it = col; it; ++it) {
            unsigned row_index = it.key() - 1; // Adjust index if needed
            matrix[row_index][i - 1] = it.val();
        }
    }

    // Print the matrix
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
}