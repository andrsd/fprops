#include "InterpolatedFluidProperties.h"
#include "Utils.h"
#include "Numerics.h"
#include "h5pp/h5pp.h"
#include "fmt/printf.h"
#include "Eigen/Dense"

namespace fprops {

const std::size_t InterpolatedFluidProperties::U = 0;
const std::size_t InterpolatedFluidProperties::V = 1;
const std::size_t InterpolatedFluidProperties::RHO = 2;
const std::size_t InterpolatedFluidProperties::P = 3;
const std::size_t InterpolatedFluidProperties::T = 4;
const std::size_t InterpolatedFluidProperties::MU = 5;
const std::size_t InterpolatedFluidProperties::CP = 6;
const std::size_t InterpolatedFluidProperties::CV = 7;
const std::size_t InterpolatedFluidProperties::S = 8;
const std::size_t InterpolatedFluidProperties::K = 9;
const std::size_t InterpolatedFluidProperties::H = 10;
const std::size_t InterpolatedFluidProperties::W = 11;

InterpolatedFluidProperties::InterpolatedFluidProperties() :
    SinglePhaseFluidProperties(),
    var_names({ "u", "v", "rho", "p", "T", "mu", "cp", "cv", "s", "k", "h", "w" }),
    var_units({ "J/kg",
                "m^3/kg",
                "kg/m^3",
                "Pa",
                "K",
                "Pa-s",
                "J/(kg-K)",
                "J/(kg-K)",
                "J/(kg-K)",
                "W/(m-K)",
                "J/kg",
                "m/s" })
{
}

void
InterpolatedFluidProperties::load(const std::string & file_name)
{
    if (std::filesystem::exists(file_name)) {
        h5pp::File file(file_name, h5pp::FileAccess::READONLY);
        this->pT_data = read_data(file, "p_T", { P, T }, { U, V, RHO, MU, CP, CV, S, K, H, W });
        this->rhoT_data = read_data(file, "rho_T", { RHO, T }, { U, V, P, MU, CP, CV, S, K, H, W });
        this->rhop_data = read_data(file, "rho_p", { RHO, P }, { T, U, V, MU, CP, CV, S, K, H, W });
        this->vu_data = read_data(file, "v_u", { V, U }, { P, T, RHO, MU, CP, CV, S, K, H, W });
        this->hs_data = read_data(file, "h_s", { H, S }, { U, V, P, T, RHO, MU, CP, CV, K, W });
    }
    else
        throw std::runtime_error(fmt::format("File '{}' does not exist.", file_name));
}

State
InterpolatedFluidProperties::p_T(double p, double T) const
{
    if (!this->pT_data.empty()) {
        auto vals = this->pT_data.eval(p, T, { RHO, U, V, MU, CP, CV, S, K, H, W });
        return { vals[U],  vals[V],  vals[RHO], p,       T,       vals[MU],
                 vals[CP], vals[CV], vals[S],   vals[K], vals[H], vals[W] };
    }
    else
        throw std::runtime_error("'p_T' data set is missing.");
}

State
InterpolatedFluidProperties::rho_T(double rho, double T) const
{
    if (!this->rhoT_data.empty()) {
        auto vals = this->rhoT_data.eval(rho, T, { P, U, V, MU, CP, CV, S, K, H, W });
        return { vals[U],  vals[V],  rho,     vals[P], T,       vals[MU],
                 vals[CP], vals[CV], vals[S], vals[K], vals[H], vals[W] };
    }
    else
        throw std::runtime_error("'rho_T' data set is missing.");
}

State
InterpolatedFluidProperties::rho_p(double rho, double p) const
{
    if (!this->rhop_data.empty()) {
        auto vals = this->rhop_data.eval(rho, p, { T, U, V, MU, CP, CV, S, K, H, W });
        return { vals[U],  vals[V],  rho,     p,       vals[T], vals[MU],
                 vals[CP], vals[CV], vals[S], vals[K], vals[H], vals[W] };
    }
    else
        throw std::runtime_error("'rho_p' data set is missing.");
}

State
InterpolatedFluidProperties::v_u(double v, double u) const
{
    if (!this->vu_data.empty()) {
        auto vals = this->vu_data.eval(v, u, { P, T, RHO, MU, CP, CV, S, K, H, W });
        return { u,        v,        vals[RHO], vals[P], vals[T], vals[MU],
                 vals[CP], vals[CV], vals[S],   vals[K], vals[H], vals[W] };
    }
    else
        throw std::runtime_error("'v_u' data set is missing.");
}

State
InterpolatedFluidProperties::h_s(double h, double s) const
{
    if (!this->hs_data.empty()) {
        auto vals = this->hs_data.eval(h, s, { U, V, P, T, RHO, MU, CP, CV, K, W });
        return { vals[U],  vals[V],  vals[RHO], vals[P], vals[T], vals[MU],
                 vals[CP], vals[CV], s,         vals[K], h,       vals[W] };
    }
    else
        throw std::runtime_error("'h_s' data set is missing.");
}

void
InterpolatedFluidProperties::check_unit(const h5pp::File & file,
                                        const std::string & dataset_name,
                                        const std::string & unit)
{
    auto attr = file.readAttribute<std::string>(dataset_name, "unit");
    if (attr != unit)
        throw std::runtime_error(fmt::format("Expected unit '{}' for '{}'.", unit, dataset_name));
}

Eigen::VectorXd
InterpolatedFluidProperties::read_var_range(const h5pp::File & file,
                                            const std::string & group_name,
                                            std::size_t var_idx)
{
    Eigen::VectorXd range;
    auto nm = fmt::format("/{}/{}", group_name, this->var_names[var_idx]);
    check_unit(file, nm, this->var_units[var_idx]);
    file.readDataset(range, nm);
    if (range.size() < 2)
        throw std::runtime_error(
            fmt::format("'{}' range must have 2 or more grid points.", this->var_names[var_idx]));
    return range;
}

DynDenseMatrix
InterpolatedFluidProperties::read_grid(const h5pp::File & file,
                                       const std::string & group_name,
                                       std::size_t var_idx)
{
    auto ds_name = fmt::format("/{}/{}", group_name, this->var_names[var_idx]);
    check_unit(file, ds_name, this->var_units[var_idx]);
    // NOTE: it should be possible to read the MatrixXd directly from the HDF5 file (there
    // are examples showing that), however I am getting compilation error on those things.
    // May be, just using wrong versions of h5cpp and Eigen.
    auto ds = file.readDataset<Eigen::VectorXd>(ds_name);
    auto dims = file.getDatasetDimensions(ds_name);
    return Eigen::Map<DynDenseMatrix>(ds.data(), dims[0], dims[1]);
}

InterpolatedFluidProperties::FPDataSet
InterpolatedFluidProperties::read_data(const h5pp::File & file,
                                       const std::string & group_name,
                                       const std::pair<std::size_t, std::size_t> & vars,
                                       const std::vector<std::size_t> & var_idx)
{
    FPDataSet fpset;
    if (file.findGroups(fmt::format("{}", group_name)).size() == 1) {
        fpset.set_range_1(read_var_range(file, group_name, vars.first));
        fpset.set_range_2(read_var_range(file, group_name, vars.second));
        for (auto & idx : var_idx) {
            auto vals = read_grid(file, group_name, idx);
            fpset.set_values(idx, vals);
        }
    }
    return fpset;
}

bool
InterpolatedFluidProperties::FPDataSet::empty() const
{
    return (this->x1.rows() == 0) || (this->x2.rows() == 0);
}

void
InterpolatedFluidProperties::FPDataSet::set_range_1(const Eigen::VectorXd & values)
{
    this->x1 = values;
}

void
InterpolatedFluidProperties::FPDataSet::set_range_2(const Eigen::VectorXd & values)
{
    this->x2 = values;
}

void
InterpolatedFluidProperties::FPDataSet::set_values(std::size_t idx, const DynDenseMatrix & y)
{
    Eigen::Matrix<double, 16, 16> weights;
    // clang-format off
    weights <<
         1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
         -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0,
         2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
         0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1,
         0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1,
         -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0,
         9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2,
         -6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2,
         2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0,
         -6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1,
         4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1;
    // clang-format on

    this->values[idx].resize(this->x1.size() - 1, this->x2.size() - 1);

    DynDenseMatrix dydx1, dydx2, d2ydx1x2;
    calc_derivatives(y, dydx1, dydx2, d2ydx1x2);

    for (Eigen::Index i = 0; i < this->values[idx].rows(); i++) {
        for (Eigen::Index j = 0; j < this->values[idx].cols(); j++) {
            const double d1 = this->x1(i + 1) - this->x1(i);
            const double d2 = this->x2(j + 1) - this->x2(j);

            Eigen::Vector4d yy0;
            yy0 << y(i, j), y(i, j + 1), y(i + 1, j + 1), y(i + 1, j);
            Eigen::Vector4d yy1;
            yy1 << dydx1(i, j), dydx1(i, j + 1), dydx1(i + 1, j + 1), dydx1(i + 1, j);
            Eigen::Vector4d yy2;
            yy2 << dydx2(i, j), dydx2(i, j + 1), dydx2(i + 1, j + 1), dydx2(i + 1, j);
            Eigen::Vector4d yy12;
            yy12 << d2ydx1x2(i, j), d2ydx1x2(i, j + 1), d2ydx1x2(i + 1, j + 1), d2ydx1x2(i + 1, j);

            Eigen::Matrix<double, 16, 1> xx;
            for (Eigen::Index k = 0; k < 4; k++) {
                xx(k) = yy0(k);
                xx(k + 4) = yy1(k) * d1;
                xx(k + 8) = yy2(k) * d2;
                xx(k + 12) = yy12(k) * d1 * d2;
            }
            Eigen::Matrix<double, 16, 1> cl = weights * xx;

            this->values[idx](i, j) = Eigen::Map<Eigen::Matrix4d>(cl.data(), 4, 4);
        }
    }
}

void
InterpolatedFluidProperties::FPDataSet::calc_derivatives(const DynDenseMatrix & y,
                                                         DynDenseMatrix & dydx1,
                                                         DynDenseMatrix & dydx2,
                                                         DynDenseMatrix & dydx1x2)
{
    const auto m = y.rows();
    const auto n = y.cols();
    dydx1.resize(m, n);
    dydx2.resize(m, n);
    dydx1x2.resize(m, n);

    for (Eigen::Index i = 0; i < m; i++) {
        for (Eigen::Index j = 0; j < n; j++) {
            // Derivative wrt x1
            if (i == 0)
                dydx1(i, j) = (y(i + 1, j) - y(i, j)) / (this->x1(i + 1) - this->x1(i));
            else if (i == m - 1)
                dydx1(i, j) = (y(i, j) - y(i - 1, j)) / (this->x1(i) - this->x1(i - 1));
            else
                dydx1(i, j) = (y(i + 1, j) - y(i - 1, j)) / (this->x1(i + 1) - this->x1(i - 1));

            // Derivative wrt x2
            if (j == 0)
                dydx2(i, j) = (y(i, j + 1) - y(i, j)) / (this->x2(j + 1) - this->x2(j));
            else if (j == n - 1)
                dydx2(i, j) = (y(i, j) - y(i, j - 1)) / (this->x2(j) - this->x2(j - 1));
            else
                dydx2(i, j) = (y(i, j + 1) - y(i, j - 1)) / (this->x2(j + 1) - this->x2(j - 1));

            // Mixed derivative d2y/dx1x2
            if (i == 0 && j == 0)
                dydx1x2(i, j) = (y(i + 1, j + 1) - y(i + 1, j) - y(i, j + 1) + y(i, j)) /
                                (this->x1(i + 1) - this->x1(i)) / (this->x2(j + 1) - this->x2(j));
            else if (i == 0 && j == n - 1)
                dydx1x2(i, j) = (y(i + 1, j) - y(i + 1, j - 1) - y(i, j) + y(i, j - 1)) /
                                (this->x1(i + 1) - this->x1(i)) / (this->x2(j) - this->x2(j - 1));
            else if (i == m - 1 && j == 0)
                dydx1x2(i, j) = (y(i, j + 1) - y(i, j) - y(i - 1, j + 1) + y(i - 1, j)) /
                                (this->x1(i) - this->x1(i - 1)) / (this->x2(j + 1) - this->x2(j));
            else if (i == m - 1 && j == n - 1)
                dydx1x2(i, j) = (y(i, j) - y(i, j - 1) - y(i - 1, j) + y(i - 1, j - 1)) /
                                (this->x1(i) - this->x1(i - 1)) / (this->x2(j) - this->x2(j - 1));
            else if (i == 0)
                dydx1x2(i, j) = (y(i + 1, j + 1) - y(i + 1, j - 1) - y(i, j + 1) + y(i, j - 1)) /
                                (this->x1(i + 1) - this->x1(i)) /
                                (this->x2(j + 1) - this->x2(j - 1));
            else if (i == m - 1)
                dydx1x2(i, j) = (y(i, j + 1) - y(i, j - 1) - y(i - 1, j + 1) + y(i - 1, j - 1)) /
                                (this->x1(i) - this->x1(i - 1)) /
                                (this->x2(j + 1) - this->x2(j - 1));
            else if (j == 0)
                dydx1x2(i, j) = (y(i + 1, j + 1) - y(i + 1, j) - y(i - 1, j + 1) + y(i - 1, j)) /
                                (this->x1(i + 1) - this->x1(i - 1)) /
                                (this->x2(j + 1) - this->x2(j));
            else if (j == n - 1)
                dydx1x2(i, j) = (y(i + 1, j) - y(i + 1, j - 1) - y(i - 1, j) + y(i - 1, j - 1)) /
                                (this->x1(i + 1) - this->x1(i - 1)) /
                                (this->x2(j) - this->x2[j - 1]);
            else
                dydx1x2(i, j) =
                    (y(i + 1, j + 1) - y(i + 1, j - 1) - y(i - 1, j + 1) + y(i - 1, j - 1)) /
                    (this->x1(i + 1) - this->x1(i - 1)) / (this->x2(j + 1) - this->x2(j - 1));
        }
    }
}

std::array<double, InterpolatedFluidProperties::FPDataSet::N>
InterpolatedFluidProperties::FPDataSet::eval(double xx1,
                                             double xx2,
                                             const std::vector<std::size_t> & val_idx) const
{
    Eigen::Index i = utils::interval_index(this->x1, xx1);
    auto u = utils::normalize_interval_location(this->x1(i), this->x1(i + 1), xx1);
    Eigen::Index j = utils::interval_index(this->x2, xx2);
    auto v = utils::normalize_interval_location(this->x2(j), this->x2(j + 1), xx2);
    std::array<double, N> vals = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for (auto & idx : val_idx) {
        const auto & coeffs = this->values[idx](i, j);
        for (Eigen::Index k = 0; k < 4; k++)
            for (Eigen::Index l = 0; l < 4; l++)
                vals[idx] += coeffs(k, l) * math::pow(u, k) * math::pow(v, l);
    }
    return vals;
}

} // namespace fprops
