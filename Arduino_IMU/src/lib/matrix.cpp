#include "matrix.h"

Matrix operator + (const double _scalar, Matrix _mat)
{
    Matrix _out(_mat.getRowNumber(), _mat.getColNumber());
    for (uint8_t _i = 0; _i < _mat.getRowNumber(); _i++)
    {
        for (uint8_t _j = 0; _j < _mat.getColNumber(); _j++)
        {
            _out[_i][_j] = _mat[_i][_j] + _scalar;
        }
    }
    return _out;
}

Matrix operator + (Matrix _mat, const double _scalar)
{
    return _scalar + _mat;
}

Matrix operator - (const double _scalar, Matrix _mat)
{
    Matrix _out(_mat.getRowNumber(), _mat.getColNumber());
    for (uint8_t _i = 0; _i < _mat.getRowNumber(); _i++)
    {
        for (uint8_t _j = 0; _j < _mat.getColNumber(); _j++)
        {
            _out[_i][_j] = -_mat[_i][_j] + _scalar;
        }
    }
    return _out;
}

Matrix operator - (Matrix _mat, const double _scalar)
{
    return -_scalar + _mat;
}

Matrix operator * (const double _scalar, Matrix _mat)
{
    Matrix _out(_mat.getRowNumber(), _mat.getColNumber());
    for (uint8_t _i = 0; _i < _mat.getRowNumber(); _i++)
    {
        for (uint8_t _j = 0; _j < _mat.getColNumber(); _j++)
        {
            _out[_i][_j] = _mat[_i][_j] * _scalar;
        }
    }
    return _out;
}

Matrix operator * (Matrix _mat, const double _scalar)
{
    return _scalar * _mat;
}
