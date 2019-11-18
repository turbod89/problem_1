#include <iostream>
#include <matrix/base.hpp>
#include <utility>
#include <set>

int look_for_solutions(int n, int argc, char * argv[]) {
    Matrix<int> m(n, n);

    // set statics
    for (int i = 1; i < n; i++) {
        m(0, i) = 1;
        m(i, 0) = 1;
    }

    //set variable
    int total_combinations = std::pow(2, (n-1) * (n-2));
    for (int signs = 0; signs < total_combinations ; signs++) {
        std::cerr << "\rMatrix " << signs + 1 << " out of " << total_combinations << std::flush;

        int c = 1;
        for (int i = 1; i < n; i++) {
            for (int j = 1; j < n; j++) {
                if (i == j) {
                    m(i,j) = 0;
                } else {
                    // get sign
                    int sign = (signs % (2*c)) / c;
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
