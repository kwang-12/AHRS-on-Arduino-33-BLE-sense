#pragma once
#include <math.h>
#include <arduino.h>
#define compare_zero (1e-7)
#define max_row_size 8
#define max_col_size 8
#define max_data_length (max_col_size * max_row_size)
/** 
 * 1. Data are stored in double (single precision)
 * 2. No boundary checks are performed. So be absolutely sure the dimensions are correct.
 **/
class Matrix
{
public:
    /**
     * basic constructor
     * fill all elements with zero
     **/
    Matrix(const int _num_row, const int _num_col)
    {
        this->_num_row = _num_row;
        this->_num_col = _num_col;
        this->setValueAll(0.0);
    }
    /**
     * constructor
     * fill matrix with given values
     **/
    Matrix(const int _num_row, const int _num_col, double *data)
    {
        this->_num_row = _num_row;
        this->_num_col = _num_col;
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                (*this)[_i][_j] = *data;
                data++;
            }
        }
    }
    int getRowNumber() { return this->_num_row; }
    int getColNumber() { return this->_num_col; }
    /**
     * make augmented matrix given the dimensions are already set 
     * as n x 2n
     * _original is a n x n matrix that needs to be transformed into an augmented matrix
     **/
    Matrix set_augment(Matrix _original)
    {
        Matrix _outp(this->_num_row, this->_num_col);
        for (int _row_tick = 0; _row_tick < this->_num_row; _row_tick++)
        {
            for (int _col_tick = 0; _col_tick < this->_num_col; _col_tick++)
            {
                if (_col_tick < this->_num_row)
                {
                    _outp[_row_tick][_col_tick] = _original[_row_tick][_col_tick];
                }
                else
                {
                    if (_row_tick == _col_tick - this->_num_row)
                    {
                        _outp[_row_tick][_col_tick] = double(1.0);
                    }
                }
            }
        }
        return _outp;
    }

    /****** operator overload *******/
    class Proxy
    {
    public:
        Proxy(double *_array, int _maxCol) : _array(_array) { this->_maxCol = _maxCol; }
        double &operator[](int _col)
        {
            return _array[_col];
        }

    private:
        double *_array;
        int _maxCol;
    };
    
    Proxy operator[](int _row)
    {
        return Proxy(_data[_row], this->_num_row);
    }
    
    bool operator == (Matrix _compare)
    {
        if ((this->_num_row != _compare._num_row) || (this->_num_col != _compare._num_col))
        {
            return false;
        }
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                if (fabs((*this)[_i][_j] - _compare[_i][_j]) > double(compare_zero))
                {
                    return false;
                }
            }
        }
        return true;
    }
    
    Matrix operator + (Matrix _add)
    {
        Matrix _out(this->_num_row, this->_num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_i][_j] = (*this)[_i][_j] + _add[_i][_j];
            }
        }
        return _out;
    }

    Matrix operator - (Matrix _sub)
    {
        Matrix _out(this->_num_row, this->_num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_i][_j] = (*this)[_i][_j] - _sub[_i][_j];
            }
        }
        return _out;
    }
    
    Matrix operator - (void)
    {
        Matrix _out(this->_num_row, this->_num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_i][_j] = -(*this)[_i][_j];
            }
        }
        return _out;
    }

    Matrix operator * (Matrix _mul)
    {
        Matrix _out(this->_num_row, _mul._num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < _mul._num_col; _j++)
            {
                _out[_i][_j] = 0.0;
                for (int _k = 0; _k < this->_num_col; _k++)
                {
                    _out[_i][_j] += ((*this)[_i][_k] * _mul[_k][_j]);
                }
            }
        }
        return _out;
    }

    Matrix operator / (const double _scalar)
    {
        Matrix _out(this->_num_row, this->_num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_i][_j] = (*this)[_i][_j] / _scalar;
            }
        }
        return _out;
    }
    
    /****** operations *******/

    /**
     * replicate matrix
     **/
    Matrix replicate()
    {
        Matrix _out(this->_num_row, this->_num_col);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_i][_j] = (*this)[_i][_j];
            }
        }
        return _out;
    }
    /** 
     * set all values
     **/
    void setValueAll(const double _val)
    {
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                (*this)[_i][_j] = _val;
            }
        }
    }
    /** 
     * set diagonal values, other values remain unchanged
     **/
    void setValueDiagonal(const double _val)
    {
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                if (_i == _j)
                {
                    (*this)[_i][_j] = _val;
                }
            }
        }
    }
    /**
     * round specific element to zero if small enough
     **/
    void roundSpecificValueToZero(const int _i, const int _j)
    {
        if (fabs((*this)[_i][_j]) < double(compare_zero))
        {
            (*this)[_j][_j] = 0.0;
        }
    }
    /**
     * round all elements to zero if small enough
     **/
    void roundMatrixToZero()
    {
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                if (fabs((*this)[_i][_j]) < double(compare_zero))
                {
                    (*this)[_j][_j] = 0.0;
                }
            }
        }
    }

    /**
     * transpose
     **/
    Matrix transpose()
    {
        Matrix _out(this->_num_col, this->_num_row);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _out[_j][_i] = (*this)[_i][_j];
            }
        }
        return _out;
    }

    /**
     * normalize vector
     * return false if the vector is a zero vector
     **/
    bool normalizeVector()
    {
        double _sum = 0.0;
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                _sum = _sum + ((*this)[_i][_j] * (*this)[_i][_j]);
            }
        }
        if (fabs(_sum) < double(compare_zero))
        {
            // is a zero vector
            return false;
        }
        _sum = sqrt(_sum);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                // normalize
                (*this)[_i][_j] /= _sum;
            }
        }
        return true;
    }

    /**
     * determinant
     * make sure the matrix is square before use
     **/
    double determinant()
    {
        double _det = 0;
        if (this->_num_row == 2) // simple 2x2
            return (((*this)[0][0] * (*this)[1][1]) - ((*this)[1][0] * (*this)[0][1]));
        else
        {
            for (int _i = 0; _i < this->_num_row; _i++)
            {
                double _submatrix_data[max_data_length];
                int _data_counter = 0;
                for (int _row_tick = 1; _row_tick < this->_num_row; _row_tick++)
                {
                    for (int _col_tick = 0; _col_tick < this->_num_col - 1; _col_tick++)
                    {
                        if (_col_tick >= _i)
                        {
                            _submatrix_data[_data_counter] = (*this)[_row_tick][_col_tick + 1];
                        }
                        else
                        {
                            _submatrix_data[_data_counter] = (*this)[_row_tick][_col_tick];
                        }
                        _data_counter++;
                    }
                }
                Matrix _submatrix(this->_num_row - 1, this->_num_col - 1, _submatrix_data);
                _det = _det + pow(-1, _i) * (*this)[0][_i] * _submatrix.determinant();
            }
        }
        return _det;
    }

    /**
     * adjoint
     **/
    Matrix adjoint()
    {
        double _out_data[max_data_length];
        int _out_counter=0;
        for (int _row_tick = 0; _row_tick < this->_num_row; _row_tick++)
        {
            for (int _col_tick = 0; _col_tick < this->_num_col; _col_tick++)
            {
                double _submatrix_data[max_data_length];
                int _data_counter = 0;
                for (int _i = 0; _i < this->_num_row-1; _i++)
                {
                    for (int _j = 0; _j < this->_num_col-1; _j++)
                    {
                        if (_i >= _row_tick)
                        {
                            if (_j >= _col_tick)
                            {
                                _submatrix_data[_data_counter] = (*this)[_i+1][_j+1];
                            }
                            else
                            {
                                _submatrix_data[_data_counter] = (*this)[_i+1][_j];
                            }
                        }
                        else
                        {
                            if (_j >= _col_tick)
                            {
                                _submatrix_data[_data_counter] = (*this)[_i][_j+1];
                            }
                            else
                            {
                                _submatrix_data[_data_counter] = (*this)[_i][_j];
                            }
                        }
                        _data_counter++;
                    }
                }
                Matrix _submatrix(this->_num_row-1, this->_num_col-1, _submatrix_data);
                _out_data[_out_counter] = pow(-1, _row_tick+_col_tick) * _submatrix.determinant();
                _out_counter++;
            }
        }
        Matrix _out(this->_num_row, this->_num_col, _out_data);
        return _out;
    }

    /**
     * invertibility check
     **/
    bool invertibilityCheck()
    {
        // check - 1
        if (this->_num_row != this->_num_col)
        {
            // non-square matrix
            return false;
        }
        // check - 2
        if (fabs(this->determinant()) < double(compare_zero))
        {
            // determinant is zero
            return false;
        }
        return true;
    }
/**
 * swap row _i with row _j
 **/
Matrix swapRow(const int _i, const int _j)
{
    Matrix _out(this->_num_row, this->_num_col);
    _out = this->replicate();
    // go through all the columns
    for (int _col_tick = 0; _col_tick < this->_num_col; _col_tick++)
    {
        double _temp = (*this)[_i][_col_tick];    // store [_ith][_col_tick] element
        // swap elements
        _out[_i][_col_tick] = _out[_j][_col_tick];
        _out[_j][_col_tick] = _temp;
    }
    return _out;
}

/**
 * inverse matrix (analytical)
 **/
Matrix inverse()
{
    Matrix _out(this->_num_row, this->_num_col);
    _out[0][0] = ((*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[1][3]*(*this)[2][2]*(*this)[3][1])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[0][1] = -((*this)[0][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][3]*(*this)[2][2]*(*this)[3][1])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[0][2] = ((*this)[0][1]*(*this)[1][2]*(*this)[3][3] - (*this)[0][1]*(*this)[1][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][1]*(*this)[3][3] + (*this)[0][2]*(*this)[1][3]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[3][2] - (*this)[0][3]*(*this)[1][2]*(*this)[3][1])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[0][3] = -((*this)[0][1]*(*this)[1][2]*(*this)[2][3] - (*this)[0][1]*(*this)[1][3]*(*this)[2][2] - (*this)[0][2]*(*this)[1][1]*(*this)[2][3] + (*this)[0][2]*(*this)[1][3]*(*this)[2][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][2] - (*this)[0][3]*(*this)[1][2]*(*this)[2][1])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);

_out[1][0] = -((*this)[1][0]*(*this)[2][2]*(*this)[3][3] - (*this)[1][0]*(*this)[2][3]*(*this)[3][2] - (*this)[1][2]*(*this)[2][0]*(*this)[3][3] + (*this)[1][2]*(*this)[2][3]*(*this)[3][0] + (*this)[1][3]*(*this)[2][0]*(*this)[3][2] - (*this)[1][3]*(*this)[2][2]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);

_out[1][1] = ((*this)[0][0]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[2][3]*(*this)[3][2] - (*this)[0][2]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[2][3]*(*this)[3][0] + (*this)[0][3]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[2][2]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[1][2] = -((*this)[0][0]*(*this)[1][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][3]*(*this)[3][2] - (*this)[0][2]*(*this)[1][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][3]*(*this)[3][0] + (*this)[0][3]*(*this)[1][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][2]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);

_out[1][3] = ((*this)[0][0]*(*this)[1][2]*(*this)[2][3] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);

_out[2][0] = ((*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[1][3]*(*this)[2][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[2][1] = -((*this)[0][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][3]*(*this)[2][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[2][2] = ((*this)[0][0]*(*this)[1][1]*(*this)[3][3] - (*this)[0][0]*(*this)[1][3]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[3][3] + (*this)[0][1]*(*this)[1][3]*(*this)[3][0] + (*this)[0][3]*(*this)[1][0]*(*this)[3][1] - (*this)[0][3]*(*this)[1][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);

_out[2][3] = -((*this)[0][0]*(*this)[1][1]*(*this)[2][3] - (*this)[0][0]*(*this)[1][3]*(*this)[2][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][3] + (*this)[0][1]*(*this)[1][3]*(*this)[2][0] + (*this)[0][3]*(*this)[1][0]*(*this)[2][1] - (*this)[0][3]*(*this)[1][1]*(*this)[2][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[3][0] = -((*this)[1][0]*(*this)[2][1]*(*this)[3][2] - (*this)[1][0]*(*this)[2][2]*(*this)[3][1] - (*this)[1][1]*(*this)[2][0]*(*this)[3][2] + (*this)[1][1]*(*this)[2][2]*(*this)[3][0] + (*this)[1][2]*(*this)[2][0]*(*this)[3][1] - (*this)[1][2]*(*this)[2][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[3][1] = ((*this)[0][0]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[2][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[3][2] = -((*this)[0][0]*(*this)[1][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[3][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
_out[3][3] = ((*this)[0][0]*(*this)[1][1]*(*this)[2][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0])/((*this)[0][0]*(*this)[1][1]*(*this)[2][2]*(*this)[3][3] - (*this)[0][0]*(*this)[1][1]*(*this)[2][3]*(*this)[3][2] - (*this)[0][0]*(*this)[1][2]*(*this)[2][1]*(*this)[3][3] + (*this)[0][0]*(*this)[1][2]*(*this)[2][3]*(*this)[3][1] + (*this)[0][0]*(*this)[1][3]*(*this)[2][1]*(*this)[3][2] - (*this)[0][0]*(*this)[1][3]*(*this)[2][2]*(*this)[3][1] - (*this)[0][1]*(*this)[1][0]*(*this)[2][2]*(*this)[3][3] + (*this)[0][1]*(*this)[1][0]*(*this)[2][3]*(*this)[3][2] + (*this)[0][1]*(*this)[1][2]*(*this)[2][0]*(*this)[3][3] - (*this)[0][1]*(*this)[1][2]*(*this)[2][3]*(*this)[3][0] - (*this)[0][1]*(*this)[1][3]*(*this)[2][0]*(*this)[3][2] + (*this)[0][1]*(*this)[1][3]*(*this)[2][2]*(*this)[3][0] + (*this)[0][2]*(*this)[1][0]*(*this)[2][1]*(*this)[3][3] - (*this)[0][2]*(*this)[1][0]*(*this)[2][3]*(*this)[3][1] - (*this)[0][2]*(*this)[1][1]*(*this)[2][0]*(*this)[3][3] + (*this)[0][2]*(*this)[1][1]*(*this)[2][3]*(*this)[3][0] + (*this)[0][2]*(*this)[1][3]*(*this)[2][0]*(*this)[3][1] - (*this)[0][2]*(*this)[1][3]*(*this)[2][1]*(*this)[3][0] - (*this)[0][3]*(*this)[1][0]*(*this)[2][1]*(*this)[3][2] + (*this)[0][3]*(*this)[1][0]*(*this)[2][2]*(*this)[3][1] + (*this)[0][3]*(*this)[1][1]*(*this)[2][0]*(*this)[3][2] - (*this)[0][3]*(*this)[1][1]*(*this)[2][2]*(*this)[3][0] - (*this)[0][3]*(*this)[1][2]*(*this)[2][0]*(*this)[3][1] + (*this)[0][3]*(*this)[1][2]*(*this)[2][1]*(*this)[3][0]);
 
    
    return _out;
}



//     /**
//      * inverse matrix (adjoint matrix method) - TOO SLOW
//      **/
//     Matrix inverse()
//     {
//         Matrix _adj = this->adjoint();
//         double _deter = this->determinant();
//         Matrix _adjT = _adj.transpose();
//         Matrix _inv = _adjT/double(_deter);
//         return _inv;
//     }
//     /**
//      * grab the inverse matrix from augmented matrix
//      **/
//     Matrix grab_inverse()
//     {
//         Matrix _out(this->_num_row, this->_num_row);
//         for (int _row_tick = 0; _row_tick < this->_num_row; _row_tick++)
//         {
//             for (int _col_tick = this->_num_row; _col_tick < this->_num_col; _col_tick++)
//             {
//                 _out[_row_tick][_col_tick-this->_num_row] = (*this)[_row_tick][_col_tick];
//             }
//         }
//         return _out;
//     }

//     // /**
//     //  * inverse matrix (gauss jordan) -  BUGGY
//     //  **/
// Matrix inverse()
// {
//     // create augmented matrix
//     Matrix _augmented(this->_num_row, this->_num_col*2);
//     _augmented = _augmented.set_augment((*this));

//     // _augmented = _augmented.InsertMatrix((*this),0,0);
//     // Matrix _eye(this->_num_row, this->_num_col);
//     // _eye.setValueDiagonal(double(1.0));
//     // _augmented = _augmented.InsertMatrix(_eye, 0, this->_num_col);
//     // interchange rows
//     for (int _i = this->_num_row-1; _i > 0; _i--)
//     {
//         if (_augmented[_i-1][0] < _augmented[_i][0])
//         {
//             // swap row i-1 with i
//             _augmented = _augmented.swapRow(_i-1, _i);
//         }
//     }
//     for (int _i = 0; _i < this->_num_row; _i++)
//     {
//         for (int _j = 0; _j < this->_num_row; _j++)
//         {
//             if (_j != _i)
//             {
//                 double _temp = _augmented[_j][_i] / _augmented[_i][_i];
//                 for (int _k = 0; _k < 2*this->_num_row; _k++)
//                 {
//                     _augmented[_j][_k] -= _augmented[_i][_k] * _temp;
//                 }
//             }
//         }
//     }
//     for (int _i = 0; _i < this->_num_row; _i++)
//     {
//         double temp = _augmented[_i][_i];
//         for (int _j = 0; _j < 2*this->_num_row; _j++)
//         {
//             _augmented[_i][_j] = _augmented[_i][_j] / temp;
//         }
//     }
//     Matrix _out = _augmented.grab_inverse();
//     return _out;
//     // _eye = _eye.InsertMatrix(_augmented,0,0,0,_augmented._num_row-1,_augmented._num_row,_augmented._num_col-1);
//     // return _eye;
// }


    /**
     * get diagonal
     * return a vertical matrix that's filled with the diagonal elements
     **/
    Matrix getDiagonal()
    {
        Matrix _diag(this->_num_row, 1);
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            _diag[_i][0] = (*this)[_i][_i];
        }
        return _diag;
    }

    /**
     * insert _vector at _pos
     * Note that the vector is vertical and its number of element should equal
     * to the array's number of rows
     * Example: A is (3,3) and B is (3,1)
     * C =A.InsertVector(B,1);
     * A = 
     * [a00 a01 a02]
     * [a10 a11 a12]
     * [a20 a21 a22]]
     * B = 
     * [b00]
     * [b10]
     * [b20]
     * C = 
     * [a00 b00 a02]
     * [a10 b10 a12]
     * [a20 b20 a22]
     **/
    Matrix InsertVector(Matrix _vector, const int _pos)
    {
        Matrix _out(this->_num_row, this->_num_col);
        _out = this->replicate();
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            _out[_i][_pos] = (*this)[_i][0];
        }
        return _out;
    }
    /**
     * insert _matrix at _posRow and _posCol
     * BE CAREFUL WITH DIMENSIONS
     * Example: A is (4,4) and B is (2,2)
     * C = A.InsertMatrix(B, 1, 1);
     * A = 
     * [a00 a01 a02 a03]
     * [a10 a11 a12 a13]
     * [a20 a21 a22 a23]
     * [a30 a31 a32 a33]]
     * B = 
     * [b00 b01]
     * [b10 b11]
     * C = 
     * [a00 a01 a02 a03]
     * [a10 b00 b01 a13]
     * [a20 b10 b11 a23]
     * [a30 a31 a32 a33]
     **/
    Matrix InsertMatrix(Matrix _matrix, const int _posRow, const int _posCol)
    {
        Matrix _out(this->_num_row, this->_num_col);
        _out = this->replicate();
        int _copy_row_len = _matrix._num_row;
        int _copy_col_len = _matrix._num_col;
        for (int _i = 0; _i < _copy_row_len; _i++)
        {
            for (int _j = 0; _j < _copy_col_len; _j++)
            {
                // Serial.print("_out[");
                // Serial.print(_i + _posRow);
                // Serial.print("][");
                // Serial.print(_j + _posCol);
                // Serial.print("]=");
                // Serial.println(_matrix[_i][_j]);
                _out[_i + _posRow][_j + _posCol] = _matrix[_i][_j];
                // Serial.println("_out is now: ");
                // _out.print();
            }
        }
        // for (int _i = 0; _i < this->_num_row; _i++)
        // {
        //     for (int _j = 0; _j < this->_num_col; _j++)
        //     {
        //         _out[_i + _posRow][_j + _posCol] = _matrix[_i][_j];
        //         _out.print();
        //     }
        // }
        return _out;
    }
    /**
     * insert a portion of _matrix at _posRow and _posCol
     * the portion is a range specified by _rowStart, _rowEnd, _colStart, _colEnd
     * Example: A is (4,4) and B is (3,3)
     * C = A.InsertMatrix(B, 1, 1, 0, 1, 0, 2)
     * A = 
     * [a00 a01 a02 a03]
     * [a10 a11 a12 a13]
     * [a20 a21 a22 a23]
     * [a30 a31 a32 a33]
     * B = 
     * [b00 b01 b02]
     * [b10 b11 b12]
     * [b20 b21 b22]
     * The section getting inserted:
     * row: 0 - 1
     * column: 0 - 2
     * note that boundary is inclusive
     * [b00 b01 b02]
     * [b10 b11 b12]
     * C =
     * [a00 a01 a02 a03]
     * [a10 b00 b01 b02]
     * [a20 b10 b11 b12]
     * [a30 a31 a32 a33]
     **/
    Matrix InsertMatrix(Matrix _matrix, const int _posRow, const int _posCol,
                        const int _rowStart, const int _rowEnd,
                        const int _colStart, const int _colEnd)
    {
        Matrix _out(this->_num_row, this->_num_col);
        _out = this->replicate();
        int _row_len = _rowEnd - _rowStart;
        int _col_len = _colEnd - _colStart;
        for (int _i = 0; _i <= _row_len; _i++)
        {
            for (int _j = 0; _j <= _col_len; _j++)
            {
                _out[_i + _posRow][_j + _posCol] = _matrix[_i + _rowStart][_j + _colStart];
            }
        }
        return _out;
    }

    void print()
    {
        char _buffer[32];
        for (int _i = 0; _i < this->_num_row; _i++)
        {
            Serial.print("[ ");
            for (int _j = 0; _j < this->_num_col; _j++)
            {
                snprintf(_buffer, sizeof(_buffer) - 1, "%e ", (*this)[_i][_j]);
                Serial.print(_buffer);
            }
            Serial.println("]");
        }
        Serial.println("");
    }

private:
    int _num_row;
    int _num_col;
    double _data[max_row_size][max_col_size] = {{0}};
};

/****** operator overload *******/
Matrix operator + (const double _scalar, Matrix _mat);
Matrix operator + (Matrix _mat, const double _scalar);
Matrix operator - (const double _scalar, Matrix _mat);
Matrix operator - (Matrix _mat, const double _scalar);
Matrix operator * (const double _scalar, Matrix _mat);
Matrix operator * (Matrix _mat, const double _scalar);
