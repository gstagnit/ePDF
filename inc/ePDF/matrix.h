//
// Authors: Rabah Abdul Khalek: rabah.khalek@gmail.com
//          Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <stdexcept>
#include <complex>

namespace ePDF
{
  /**
   * @name Class to manage matrix objects
   */
  ///@{
  template <class T>
  class Matrix
  {
  public:
    //_________________________________________________________________________________
    Matrix(int const& Lines = 0, int const& Columns = 0, std::vector<T> const& Entries = {}):
      _Lines(Lines),
      _Columns(Columns),
      _Matrix(Entries)
    {
      // Check that the size of the Entries match the size of the matrix
      if (_Lines * _Columns != (int) Entries.size())
        throw std::runtime_error("[Matrix::Matrix]: the size of the input vector does not match the size of the matrix.");
    }

    //_________________________________________________________________________________
    void SetElement(int const& i, int const& j, T const& value)
    {
      if (i < 0 || i > _Lines)
        throw std::runtime_error("[Matrix::SetElement]: line index out of range.");

      if (j < 0 || j > _Columns)
        throw std::runtime_error("[Matrix::SetElement]: column index out of range.");

      _Matrix[i * _Columns + j] = value;
    }

    //_________________________________________________________________________________
    T GetElement(int const& i, int const& j) const
    {
      if (i < 0 || i > _Lines)
        throw std::runtime_error("[Matrix::GetElement]: line index out of range.");

      if (j < 0 || j > _Columns)
        throw std::runtime_error("[Matrix::GetElement]: column index out of range.");

      return _Matrix[i * _Columns + j];
    }

    //_________________________________________________________________________________
    T operator () (int const& i, int const& j) const
    {
      return GetElement(i, j);
    }

    //_________________________________________________________________________________
    int GetLines() const
    {
      return _Lines;
    }

    //_________________________________________________________________________________
    int GetColumns() const
    {
      return _Columns;
    }

    //_________________________________________________________________________________
    std::vector<T> GetVector() const
    {
      return _Matrix;
    }

    //_________________________________________________________________________________
    Matrix<T> operator = (Matrix<T> const& term)
    {
      _Lines = term.GetLines();
      _Columns = term.GetColumns();
      _Matrix = term.GetVector();
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator += (Matrix<T> const& term)
    {
      if (_Lines != term.GetLines() || _Columns != term.GetColumns())
        throw std::runtime_error("[Matrix::operator +=]: Lines or Columns don't match adding the two matrices.");

      const std::vector<T> v = term.GetVector();
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] += v[i * _Columns + j];
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator -= (Matrix<T> const& term)
    {
      if (_Lines != term.GetLines() || _Columns != term.GetColumns())
        throw std::runtime_error("[Matrix::operator -=]: Lines or Columns don't match adding the two matrices.");

      const std::vector<T> v = term.GetVector();
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] -= v[i * _Columns + j];
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator *= (Matrix<T> const& term)
    {
      const int l1 = _Lines;
      const int c1 = _Columns;
      const int l2 = term.GetLines();
      const int c2 = term.GetColumns();
      if (c1 != l2)
        throw std::runtime_error("[Matrix::operator *=]: Lines or Columns don't match multiplying the two matrices.");

      const std::vector<T> v = term.GetVector();
      std::vector<T> out(l1 * c2, 0.);
      for (int i = 0; i < l1; i++)
        for (int j = 0; j < c2; j++)
          for (int k = 0; k < c1; k++)
            out[i * c2 + j] += _Matrix[i * c1 + k] * v[k * c2 + j];

      _Lines = l1;
      _Columns = c2;
      _Matrix = out;
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator *= (T const& coef)
    {
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] *= coef;
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator /= (T const& coef)
    {
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] /= coef;
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator *= (double const& coef)
    {
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] *= coef;
      return *this;
    }

    //_________________________________________________________________________________
    Matrix<T> operator /= (double const& coef)
    {
      for (int i = 0; i < _Lines; i++)
        for (int j = 0; j < _Columns; j++)
          _Matrix[i * _Columns + j] /= coef;
      return *this;
    }

  private:
    int _Lines;
    int _Columns;
    std::vector<T> _Matrix;

    template<class Y>
    friend std::ostream& operator << (std::ostream& os, Matrix<Y> const& sg);
  };

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, Matrix<std::complex<double>> const& m);

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator + (Matrix<T> lhs, Matrix<T> const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator - (Matrix<T> lhs, Matrix<T> const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator * (Matrix<T> lhs, Matrix<T> const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator * (Matrix<T> lhs, T const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator * (T const& lhs, Matrix<T> rhs)
  {
    return rhs *= lhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator / (Matrix<T> lhs, T const& rhs)
  {
    return lhs /= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator / (T const& lhs, Matrix<T> rhs)
  {
    return rhs /= lhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator * (Matrix<T> lhs, double const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator * (double const& lhs, Matrix<T> rhs)
  {
    return rhs *= lhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator / (Matrix<T> lhs, double const& rhs)
  {
    return lhs /= rhs;
  }

  //_________________________________________________________________________
  template<class T>
  Matrix<T> operator / (double const& lhs, Matrix<T> rhs)
  {
    return rhs /= lhs;
  }
  ///@}
}
