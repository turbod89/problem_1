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

        Matrix<int> k = m.kernel();
        m.transpose();
        std::pair<Matrix<int>, Matrix<int>> p = m.gauss_seidel();

        if ( !check_matrix_is_triangular(p.first) ) {
            std::cout << "Matrix" << std::endl;
            std::cout << m << std::endl;
            std::cout << "has Triangualation" << std::endl;
            std::cout << p.first << std::endl;
            std::cout << "has Transformation" << std::endl;
            std::cout << p.second << std::endl;
        }
        m.transpose();

        if (k.size(1) > 0) {
            k.transpose();

            if (k.size(0) > 1) {
                std::cout << "This kernel has dimension " << k.size(0) << std::endl;
                std::cout << "Matrix" << std::endl;
                std::cout << m << std::endl;
                std::cout << "has kernel" << std::endl;
                std::cout << k << std::endl;
                return 0;
            }

            for (int i = 0; i < k.size(0); i++) {
                // check all components are different in abs value
                std::set<int> components;
                bool already_in = false;
                for (int j = 0; j < k.size(1) && !already_in; j++) {
                    already_in = components.find(abs(k(i,j))) == components.end();
                }

                if (!already_in) {
                    std::cout << "Matrix" << std::endl;
                    std::cout << m << std::endl;
                    std::cout << "has kernel" << std::endl;
                    std::cout << k << std::endl;
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

    /*
    Matrix<int> m;
    std::cin >> m;

    std::cout << "Matrix" << std::endl;
    std::cout << m << std::endl;
    m.transpose();
    std::pair<Matrix<int>, Matrix<int> > p = m.gauss_seidel();
    m.transpose();
    Matrix<int> k = m.kernel();
    std::cout << "has Triangualation" << std::endl;
    std::cout << p.first << std::endl;
    std::cout << "has Transformation" << std::endl;
    std::cout << p.second << std::endl;
    std::cout << "has kernel" << std::endl;
    std::cout << k << std::endl;
    */

    return 0;    
}
