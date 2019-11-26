#include <iostream>
#include <matrix/base.hpp>
#include <utility>
#include <set>

std::string time_interval_to_string(const double & t) {
    return
        std::to_string((int)(t/3600)) + "h "
        + std::to_string((int)(t/60) % 60) + "m "
        + std::to_string((int)(t/60) % 60) + "s"
    ;
}

int _gcd(int * coefs, int n) {
    int min = 0, max = 0;
    int index_min = -1;

    // take the min (!= 0) and max
    for (int i = 0; i < n; i++) {
        if (coefs[i] != 0 ) {
            if (index_min == -1 || min > coefs[i]) {
                min = coefs[i];
                index_min = i;
            }

            if (max < coefs[i]) {
                max = coefs[i];
            }
        }
    }

    if (index_min < 0) {
        return 0;
    }

    if (max == min) {
        return max;
    }

    if (min == 1) {
        return 1;
    }

    for (int i = 0; i < n; i++) {
        if (i != index_min) {
            coefs[i] = coefs[i] % min;
        }
    }

    return _gcd(coefs, n);
}

int gcd(const int * coefs, int n) {
    int * c = (int *) std::malloc(sizeof(int) * n);
    // std::memcpy(c, coefs, n);
    for (int i = 0; i < n; i++) {
        c[i] = std::abs(coefs[i]);
    }
    int g = _gcd(c, n);
    std::free(c);
    return g;
}

bool check_matrix_is_triangular(const Matrix<int> & m) {
    int n = std::min(m.size(0), m.size(1));
    for (int i = 0; i < m.size(0); i++ ) {
        for (int j = 0; j < std::min(n, i); j++) {
            if (m(i,j) != 0) {
                return false;
            }
        }
    }
    return true;
}

void change_column_sign(Matrix<int>& m, int j) {
    for (int i = 0; i < m.size(0); i++) {
        m(i, j) *= -1;
    }
}

bool check_coefs_are_non_zero_and_different(const Matrix<int> & m, int i) {

    // check all components are different in abs value
    std::set<int> components;
    bool already_in = false;
    for (int j = 0; j < m.size(1) && !already_in; j++) {
        if (m(i, j) == 0) {
            return false;
        }
        already_in = components.find(std::abs(m(i,j))) != components.end();
        components.insert(abs(m(i, j)));
    }

    return !already_in;
}

int look_for_solutions(int n, int argc, char * argv[]) {

    time_t start_time, end_time, last_stats_update_time;
    double time_taken = 0, estimated_time = 0, time_since_last_stats_update = 0;
    const double stats_update_interval = 30;
    time(&start_time);
    time(&last_stats_update_time);

    Matrix<int> m(n, n);
    
    // set statics
    for (int i = 1; i < n; i++) {
        m(0, i) = 1;
        m(i, 0) = 1;
    }

    //set variable
    long int total_combinations = std::pow(2, (n-1) * (n-2));
    std::cerr << "Calculing kernels of " << total_combinations << " matrices." << std::endl;
    for (long int signs = 0; signs < total_combinations ; signs++) {

        time(&end_time);
        time_taken = double(end_time - start_time);
        time_since_last_stats_update = double(end_time - last_stats_update_time);

        if (time_since_last_stats_update >= stats_update_interval ) {

            estimated_time = time_taken * double(total_combinations - signs) / double(signs);

            std::cerr << "\rMatrix " << signs + 1 << " out of " << total_combinations << ".";
            std::cerr.precision(4);
            std::cerr << " ( " << std::fixed << 100 * double(signs) / double(total_combinations) << "% ).";
            std::cerr << " Estimated time to finish: ";
            std::cerr << " " << time_interval_to_string(estimated_time) << ".";
            std::cerr << std::flush;
            time(&last_stats_update_time);
        }

        // Build all -1 / +1
        long int c = 1;
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                if (i == j) {
                    m(i,j) = 0;
                } else {
                    // get sign
                    long int sign = (signs % (2*c)) / c;
                    if (sign == 0) {
                        m(i, j) = -1;
                    } else {
                        m(i, j) = 1;
                    }
                    c *= 2;
                }
            }
        }

        // calculate kernel
        Matrix<int> k = m.kernel();

        // check kernel
        if (k.size(1) > 0) {
            k.transpose();

            for (int i = 0; i < k.size(0); i++) {
                // NOTE: if k.size(0) > 1, we could be skiping possible results
            
                // Check kernel
                if (check_coefs_are_non_zero_and_different(k, i)) {

                    // copy coefs
                    int * coefs = (int *) std::malloc(sizeof(int) * n);
                    std::memcpy(coefs, &(k(i,0)), sizeof(int) * n);

                    // turn them positive (and adjust matrix consequently)
                    for (int j = 0; j < n; j++) {
                        if (coefs[j] < 0) {
                            change_column_sign(m, j);
                            coefs[j] *= -1;
                        }
                    }

                    // divide by gcd (to get minimal solution)
                    int g = gcd(coefs, n);
                    for(int j = 0 ; j < n; j++) {
                        coefs[j] /= g;
                    }

                    // get coefs total sum
                    int s = 0;
                    for (int j = 0 ; j < n ; j++) {
                        s += coefs[j];
                    }

                    // print coefficients
                    std::cout << std::endl;
                    std::cout << "{ ";
                    for(int j = 0 ; j < n; j++) {
                        if (j > 0) {
                            std::cout << ", ";
                         }
                         std::cout << coefs[j];
                    }
                    std::cout << " }" << std::endl;

                    // print equations
                    for (int l = 0; l < m.size(0); l++) {

                        std::cout << (s - coefs[l])/2 << " = ";

                        bool is_first = true;
                        for (int j = 0; j < m.size(1); j++) {
                            if (m(l,j) > 0) {
                                if (!is_first) {
                                    std::cout << " + ";
                                }
                                std::cout << coefs[j];
                                is_first = false;
                            }
                        }

                        is_first = true;
                        std::cout << " = ";
                        for (int j = 0; j < m.size(1); j++) {
                            if (m(l,j) < 0) {
                                if (!is_first) {
                                    std::cout << " + ";
                                }
                                std::cout << coefs[j];
                                is_first = false;
                            }
                        }
                        std::cout << std::endl;
                    }

                    // free memory
                    std::free(coefs);
                }
            }

        }

    }

    time(&end_time);
    time_taken = double(end_time - start_time);
    std::cerr << "It took a total time of " << time_interval_to_string(time_taken) << ".";
    std::cerr << std::endl;

    return -1;
}

int main (int argc, char * argv[]) {

    for (int n = 3; n < 9; n += 2) {
        std::cerr << "Testing for n = " << n << std::endl;
        int status = look_for_solutions(n, argc, argv);
        if (status == 0) {
            return status;
        }
    }

    return 0;
}
