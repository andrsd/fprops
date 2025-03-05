// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "fprops/single_phase_fluid_properties.h"
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <string>

namespace h5pp {
class File;
}

namespace fprops {

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DynDenseMatrix;
typedef Eigen::Matrix<Eigen::Matrix4d, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    CoefficientMatrix;

class InterpolatedFluidProperties : public SinglePhaseFluidProperties {
public:
    InterpolatedFluidProperties();
    virtual ~InterpolatedFluidProperties() = default;

    /// Load data set from an HDF5 file
    ///
    /// @param file_name HDF5 file name
    virtual void load(const std::string & file_name);

    State p_T(double p, double T) const;
    State rho_T(double rho, double T) const;
    State rho_p(double rho, double p) const;
    State v_u(double v, double u) const;
    State h_s(double h, double s) const;

protected:
    /// Fluid property data set
    ///
    /// Data set consists of 2 independent variables and a set of dependent values.
    /// Independent variables are given as an 1D array of values in ascending order.
    /// Dependent variables are given as a 2D array of values with dimension (size of 1st indep.
    /// variable, size od 2nd indep. variable).
    class FPDataSet {
        /// Number of dependent values
        static constexpr std::size_t N = 12;

    public:
        /// Check if the data exists
        ///
        /// @return `true` if this is a valid data set, `false` is this is an empty data set
        bool empty() const;

        /// Evaluate variables at location (x1, x2)
        ///
        /// @param x1 Value of 1st independent variable
        /// @param x2 Value of 2nd independent variable
        /// @param val_idx Indices of dependent variables to compute
        /// @return Values of computed variables
        std::array<double, N>
        eval(double x1, double x2, const std::vector<std::size_t> & val_idx) const;

        /// Set range for 1st independent variable
        ///
        /// @param values Values of 1st independent variable
        void set_range_1(const Eigen::VectorXd & values);

        /// Set range for 2nd independent variable
        ///
        /// @param values Values of 2nd independent variable
        void set_range_2(const Eigen::VectorXd & values);

        /// Compute coefficient matrix for a dependent variable
        ///
        /// @param idx Index of the dependent variable
        /// @param y Data grid with dependent variable values
        void set_values(std::size_t idx, const DynDenseMatrix & y);

    protected:
        /// List of the first independent variable
        Eigen::VectorXd x1;
        /// List of the second independent variable
        Eigen::VectorXd x2;
        /// Coefficient matrices for dependent variables
        std::array<CoefficientMatrix, N> values;

        /// Calculate derivatives of a 2D array with variable values
        ///
        /// @param y Variable values
        /// @param dydx1 Computed 1st derivative wrt 1st independent variable
        /// @param dydx2 Computed 1st derivative wrt 2st independent variable
        /// @param dydx1x2 Computed mixed derivative wrt 1st and 2nd independent variable
        void calc_derivatives(const DynDenseMatrix & y,
                              DynDenseMatrix & dydx1,
                              DynDenseMatrix & dydx2,
                              DynDenseMatrix & dydx1x2);
    };

    /// Check units listed on the dataset
    ///
    /// If expected unit is not found and exception is thrown
    ///
    /// @param file HDF5 file
    /// @param dataset_name Dataset name
    /// @param unit Expected unit
    void
    check_unit(const h5pp::File & file, const std::string & dataset_name, const std::string & unit);

    /// Read fluid property values from a group
    ///
    /// @param file HDF5 file
    /// @param group_name Group name to read from
    /// @param vars Independent variables
    /// @param idxs Indices of dependent variables
    /// @return Fluid property data set to evaluate as a list of variable indices
    FPDataSet read_data(const h5pp::File & file,
                        const std::string & group_name,
                        const std::pair<std::size_t, std::size_t> & vars,
                        const std::vector<std::size_t> & idxs);

    /// Read independent variable values
    ///
    /// @param file HDF5 file
    /// @param group_name Group name
    /// @param var_idx Index of an independent variable
    /// @return Vector of independent variable values
    Eigen::VectorXd
    read_var_range(const h5pp::File & file, const std::string & group_name, std::size_t var_idx);

    /// Read dependent variable values
    ///
    /// @param file HDF5 file
    /// @param group_name Group name
    /// @param var_idx Index of an dependent variable
    /// @return Grid values as a 2D array (represented as a matrix)
    DynDenseMatrix
    read_grid(const h5pp::File & file, const std::string & group_name, std::size_t var_idx);

    /// Indices for computed quantities and indexing into var_names, var_units, etc.
    static const std::size_t U;
    static const std::size_t V;
    static const std::size_t RHO;
    static const std::size_t P;
    static const std::size_t T;
    static const std::size_t MU;
    static const std::size_t CP;
    static const std::size_t CV;
    static const std::size_t S;
    static const std::size_t K;
    static const std::size_t H;
    static const std::size_t W;

    /// Fluid properties depending on pressure and temperature
    FPDataSet pT_data;
    /// Fluid properties depending on density and temperature
    FPDataSet rhoT_data;
    /// Fluid properties depending on density and pressure
    FPDataSet rhop_data;
    /// Fluid properties depending on specific volume and internal energy
    FPDataSet vu_data;
    /// Fluid properties depending on enthalpy and entropy
    FPDataSet hs_data;
    /// Variable names
    std::vector<std::string> var_names;
    /// Variable units
    std::vector<std::string> var_units;
};

} // namespace fprops
