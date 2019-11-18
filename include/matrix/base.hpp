#ifndef INCLUDE_MATRIX_
#define INCLUDE_MATRIX_

#include <iostream>
#include <cstring>
#include <sstream>
#include <cassert>
#include <cmath>
#include <utility>

template<class Ring>
class Matrix {

    private:
        int * _size;
        Ring * _values;

    public:

        ~Matrix() {
            std::free(this->_size);
            std::free(this->_values);
        }

        Matrix() {
            this->_size = new (std::nothrow) int [2];
            this->_size[0] = 0;
            this->_size[1] = 0;
            this->_values = nullptr;
        }

        Matrix(int n, int m) {
            this->_size = new (std::nothrow) int [2];
            this->_size[0] = n;
            this->_size[1] = m;
            this->_values = new (std::nothrow) Ring [this->size()];
            
            std::memset(this->_values, 0, sizeof(Ring) * this->size());
        }

        Matrix(const Matrix & M) {
            this->_size = new (std::nothrow) int [2];
            this->_size[0] = M.size(0);
            this->_size[1] = M.size(1);
            this->_values = new (std::nothrow) Ring [this->size()];
            std::memcpy(this->_values, M._values, sizeof(Ring) * this->size());
        }

        virtual int size() const {
            return this->_size[0] * this->_size[1];
        }

        virtual int size(int i) const {
            return this->_size[i];
        }

        virtual const Ring & get(int i, int j) const {
            return this->_values[i * this->size(1) + j];
        }

        virtual Ring & get(int i, int j) {
            return const_cast<int &>(static_cast<const Matrix &>(*this).get(i, j));
        }

        static Matrix Identity(int n) {
            Matrix m = Matrix(n, n);
            for (int i = 0; i < n; i++) {
                m(i, i) = (Ring) 1;
            }
            return m;
        }

        // Read & Write
        std::ostream & write(std::ostream & out) const {
            //out << "[ ";
            out << this->size(0) << " " << this->size(1) << std::endl;
            for (int i = 0 ; i < this->size(0); i++) {

                if (i > 0) {
                    out << std::endl;
                    // out << "  ";
                }

                
                for (int j = 0; j < this->size(1); j++) {
                    if (j > 0) {
                        out << " ";
                    }
                    
                    out << this->get(i, j);
                }
            }
            //out << " ]";
            out << std::endl;
        };

        friend std::ostream & operator<<(std::ostream& out, const Matrix & m) {
            return m.write(out);
        }

        std::istream & read(std::istream & in) {
            in >> this->_size[0] >> this->_size[1];
            this->_values = (Ring *) std::realloc(this->_values, sizeof(Ring) * this->size());

            for (int k = 0; k < this->size(); k++) {
                in >> this->_values[k];
            }

            return in;
        }

        friend std::istream & operator>>(std::istream & in, Matrix & m) {
            return m.read(in);
        }

        // Overloads
        const Ring & operator()(int i, int j) const {
            return this->get(i, j);
        }

        Ring & operator()(int i, int j) {
            return const_cast<Ring &>(static_cast<const Matrix &>(*this).operator()(i, j));
        }

        Matrix operator+(const Matrix & a) const {
            return Matrix::add( * this, a);
        }

        // Matrix operations
        Matrix & swap_rows(int i, int j) {
            void * aux_memory = std::malloc(sizeof(Ring) * this->size(1));
            std::memcpy(
                aux_memory,
                (void *) &(this->_values[i*this->size(1)]),
                sizeof(Ring) * this->size(1)
            );
            std::memcpy(
                (void *) &(this->_values[i*this->size(1)]),
                (void *) &(this->_values[j*this->size(1)]),
                sizeof(Ring) * this->size(1)
            );
            std::memcpy(
                (void *) &(this->_values[j*this->size(1)]),
                aux_memory,
                sizeof(Ring) * this->size(1)
            );
            std::free(aux_memory);
            return (*this);
        }

        Matrix & swap_cols(int i, int j) {
            Ring aux;
            for (int k = 0; k < this->size(0); k++) {
                aux = this->get(k, i);
                this->get(k, i) = this->get(k, j);
                this->get(k, j) = aux;
            }
            return (* this);
        }

        Matrix & add(const Matrix & m) {
            assert(m.size(0) == this->size(0));
            assert(m.size(1) == this->size(1));
            for (int k = 0; k < this->size(); k++) {
                this->_values[k] += m._values[k];
            }

            return (* this);
        }

        static Matrix add(const Matrix & a, const Matrix & b) {
            Matrix c(a);
            return c.add(b);
        }

        Matrix & transpose() {
            // NOTE: Of course, it can be optimized

            Ring * values = (Ring *) std::malloc(sizeof(Ring) * this->size());
            for (int i = 0; i < this->size(0) ; i ++) {
                for (int j = 0 ; j < this->size(1) ; j ++) {
                    values[j * this->size(0) + i] = this->_values[i * this->size(1) + j];
                }
            }
            std::memcpy(this->_values, values, sizeof(Ring) * this->size());
            std::free(values);

            int aux = this->_size[0];
            this->_size[0] = this->_size[1];
            this->_size[1] = aux;

            return (* this);
        }

        static Matrix transpose(const Matrix & a) {
            Matrix c(a);
            return c.transpose();
        }

        // Gauss-Seidel
        std::pair<Matrix, Matrix> gauss_seidel() const {
            // NOTE: It can be optimized by looking for the
            // GCD of the pivots

            int n = std::min(this->size(0), this->size(1));
            Matrix a(*this);
            Matrix b = Matrix::Identity(this->size(0));
            

            for (int i = 0; i < n; i++) {
                // get row with minor, non-zero, coefficient

                int selected_row_index = i;
                Ring selected_row_pivot_value = a(i,i);
                for (int k = i+1; k < a.size(0); k++) {
                    if (
                        a.get(k, i) != (Ring) 0
                        && (
                            selected_row_pivot_value == (Ring) 0
                            || std::abs(selected_row_pivot_value) > std::abs(a.get(k, i))
                        )
                    ) {
                        selected_row_index = k;
                        selected_row_pivot_value = a.get(k, i);
                    }
                }

                if (selected_row_index != i ) {
                    a.swap_rows(i, selected_row_index);
                    b.swap_rows(i, selected_row_index);
                }


                // for each row, apply the cross operation
                for (int k = i+1; k < a.size(0) ; k++) {
                    // apply
                    Ring row_pivot = a(k, i);
                    for (int j = i; j < a.size(1) ; j++) {
                        // NOTE: this can be optimized by starting by j = i + 1 
                        // and set the whole column at 0.

                        a(k, j) = a(i,i) * a(k, j) - row_pivot * a(i, j);
                    }
                    // apply to transformation matrix
                    for (int j = 0; j < a.size(0) ; j++) {
                        // NOTE: this can be optimized by starting by j = i + 1 
                        // and set the whole column at 0.

                        b(k, j) = a(i,i) * b(k, j) - row_pivot * b(i, j);
                    }
                }
            }


            std::pair<Matrix, Matrix> p(a, b);
            return p;
        }

        // Get kernel basis
        Matrix kernel() const {
            Matrix c(* this);
            c.transpose();
            std::pair<Matrix, Matrix > p = c.gauss_seidel();

            // loop over rows searching for a zero row
            int k = 0;
            bool all_zeros = false;
            while ( k < p.first.size(0) && !all_zeros) {
                all_zeros = true;
                for (int j = 0; j < p.first.size(1) && all_zeros; j++) {
                    all_zeros = all_zeros && p.first(k, j) == (Ring) 0;
                }
                k++;
            }
            k--;
            
            if (!all_zeros) {
                Matrix a(0,0);
                return a;
            }

            
            p.second.transpose();
            // select from k (included) until last column of p.second.transpose()

            Matrix a(p.second.size(0) , p.second.size(1) - k);
            for (int i = 0; i < a.size(0); i++) {
                for (int j = 0; j < a.size(1); j++) {
                    a(i, j) = p.second(i, j + k);
                }
            }
            
            return a;
        }

};

#endif