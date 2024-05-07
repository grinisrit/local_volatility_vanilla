/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */

#pragma once

#include <vector>


// spline interpolation
class spline
{

public:

    // default constructor: set boundary condition to be zero curvature
    // at both ends, i.e. natural splines
    spline(): 
        m_c0(0.0), m_left_value(0.0)
    {
        ;
    }
   
    spline(const std::vector<double>& X, const std::vector<double>& Y)
    {
        //m_type = cspline;
        //m_left = second_deriv;
        //m_right = second_deriv;
        //m_left_value = 0;
        //m_right_value = 0;
        //this->set_points(X,Y);
        const size_t n = X.size();
        m_x.reserve(n);
        for(const auto& x: X) m_x.push_back(x);
        m_y.reserve(n);
        for(const auto& y: Y) m_y.push_back(y);

        m_b.reserve(n);
        for(const auto& x: X) m_b.push_back(x);
        m_c.reserve(n);
        for(const auto& x: X) m_c.push_back(x);
        m_d.reserve(n);
        for(const auto& x: X) m_d.push_back(x);
    }


    // set all data points (cubic_spline=false means linear interpolation)
    
    //void set_points(const std::vector<double>& x, const std::vector<double>& y);


    // evaluates the spline at point x
    //double operator() (double x) const;

    //double get_x_min() const { assert(!m_x.empty()); return m_x.front(); }
    //double get_x_max() const { assert(!m_x.empty()); return m_x.back(); }

private:
    double m_c0;                            // for left extrapolation
    double m_left_value;
    std::vector<double> m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i, or else it won't go through grid points
    std::vector<double> m_b,m_c,m_d;        // spline coefficients

};


/*
namespace internal
{

// band matrix solver
class band_matrix
{
private:
    std::vector< std::vector<double> > m_upper;  // upper band
    std::vector< std::vector<double> > m_lower;  // lower band
public:
    band_matrix() {};                             // constructor
    band_matrix(int dim, int n_u, int n_l);       // constructor
    ~band_matrix() {};                            // destructor
    void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension
    int num_upper() const
    {
        return (int)m_upper.size()-1;
    }
    int num_lower() const
    {
        return (int)m_lower.size()-1;
    }
    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diagonal (in m_lower)
    double& saved_diag(int i);
    double  saved_diag(int i) const;
    void lu_decompose();
    std::vector<double> r_solve(const std::vector<double>& b) const;
    std::vector<double> l_solve(const std::vector<double>& b) const;
    std::vector<double> lu_solve(const std::vector<double>& b,
                                 bool is_lu_decomposed=false);

};

} // namespace internal
*/



// ---------------------------------------------------------------------
// implementation part, which could be separated into a cpp file
// ---------------------------------------------------------------------

// spline implementation
// -----------------------




/*
namespace internal
{

// band_matrix implementation
// -------------------------

band_matrix::band_matrix(int dim, int n_u, int n_l)
{
    resize(dim, n_u, n_l);
}
void band_matrix::resize(int dim, int n_u, int n_l)
{
    assert(dim>0);
    assert(n_u>=0);
    assert(n_l>=0);
    m_upper.resize(n_u+1);
    m_lower.resize(n_l+1);
    for(size_t i=0; i<m_upper.size(); i++) {
        m_upper[i].resize(dim);
    }
    for(size_t i=0; i<m_lower.size(); i++) {
        m_lower[i].resize(dim);
    }
}
int band_matrix::dim() const
{
    if(m_upper.size()>0) {
        return m_upper[0].size();
    } else {
        return 0;
    }
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
    if(k>=0)    return m_upper[k][i];
    else        return m_lower[-k][i];
}
double band_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
    if(k>=0)    return m_upper[k][i];
    else        return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double band_matrix::saved_diag(int i) const
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}
double & band_matrix::saved_diag(int i)
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose()
{
    int  i_max,j_max;
    int  j_min;
    double x;

    // preconditioning
    // normalize column i so that a_ii=1
    for(int i=0; i<this->dim(); i++) {
        assert(this->operator()(i,i)!=0.0);
        this->saved_diag(i)=1.0/this->operator()(i,i);
        j_min=std::max(0,i-this->num_lower());
        j_max=std::min(this->dim()-1,i+this->num_upper());
        for(int j=j_min; j<=j_max; j++) {
            this->operator()(i,j) *= this->saved_diag(i);
        }
        this->operator()(i,i)=1.0;          // prevents rounding errors
    }

    // Gauss LR-Decomposition
    for(int k=0; k<this->dim(); k++) {
        i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
        for(int i=k+1; i<=i_max; i++) {
            assert(this->operator()(k,k)!=0.0);
            x=-this->operator()(i,k)/this->operator()(k,k);
            this->operator()(i,k)=-x;                         // assembly part of L
            j_max=std::min(this->dim()-1,k+this->num_upper());
            for(int j=k+1; j<=j_max; j++) {
                // assembly part of R
                this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
            }
        }
    }
}
// solves Ly=b
std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_start;
    double sum;
    for(int i=0; i<this->dim(); i++) {
        sum=0;
        j_start=std::max(0,i-this->num_lower());
        for(int j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
        x[i]=(b[i]*this->saved_diag(i)) - sum;
    }
    return x;
}
// solves Rx=y
std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_stop;
    double sum;
    for(int i=this->dim()-1; i>=0; i--) {
        sum=0;
        j_stop=std::min(this->dim()-1,i+this->num_upper());
        for(int j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
        x[i]=( b[i] - sum ) / this->operator()(i,i);
    }
    return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
        bool is_lu_decomposed)
{
    assert( this->dim()==(int)b.size() );
    std::vector<double>  x,y;
    if(is_lu_decomposed==false) {
        this->lu_decompose();
    }
    y=this->l_solve(b);
    x=this->r_solve(y);
    return x;
}

} // namespace internal
*/
