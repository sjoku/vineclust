#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
double IFM_margin_obj_cpp(NumericVector pars, NumericVector data_p, NumericVector z_values, std::string family) {
    double par1 = 0, par2 = 0, par3 = 0, par4 = 0;
    
    // Regularize parameters based on family
    if (family == "Normal" || family == "Lognormal" || family == "Logistic" || family == "Cauchy") {
        par1 = pars[0];
        par2 = std::max(1e-10, pars[1]);
    } else if (family == "Skew Normal" || family == "Student-t") {
        par1 = pars[0];
        par2 = std::max(1e-10, pars[1]);
        par3 = pars[2];
    } else if (family == "Skew Student-t") {
        par1 = pars[0];
        par2 = std::max(1e-10, pars[1]);
        par3 = pars[2];
        par4 = pars[3];
    } else {
        // Gamma and Loglogistic: both parameters strictly positive
        par1 = std::max(1e-10, pars[0]);
        par2 = std::max(1e-10, pars[1]);
    }
    
    int n = data_p.size();
    NumericVector m_dens(n);
    
    if (family == "Cauchy") {
        for(int i = 0; i < n; i++) m_dens[i] = R::dcauchy(data_p[i], par1, par2, 0);
    } else if (family == "Gamma") {
        for(int i = 0; i < n; i++) m_dens[i] = R::dgamma(data_p[i], par1, 1.0/par2, 0);
    } else if (family == "Logistic") {
        for(int i = 0; i < n; i++) m_dens[i] = R::dlogis(data_p[i], par1, par2, 0);
    } else if (family == "Loglogistic") {
        for(int i = 0; i < n; i++) {
            if (data_p[i] > 0) {
                double term = std::pow(data_p[i] * par2, par1);
                m_dens[i] = (par1 * par2 * std::pow(data_p[i] * par2, par1 - 1.0)) / std::pow(1.0 + term, 2.0);
            } else {
                m_dens[i] = 0.0;
            }
        }
    } else if (family == "Lognormal") {
        for(int i = 0; i < n; i++) m_dens[i] = R::dlnorm(data_p[i], par1, par2, 0);
    } else if (family == "Normal") {
        for(int i = 0; i < n; i++) m_dens[i] = R::dnorm(data_p[i], par1, par2, 0);
    } else if (family == "Skew Normal" || family == "Student-t" || family == "Skew Student-t") {
        // Fallback to R function using Rcpp::Function
        Environment fGarch = Environment::namespace_env("fGarch");
        if (family == "Skew Normal") {
            Function dsnorm = fGarch["dsnorm"];
            m_dens = dsnorm(data_p, _["mean"] = par1, _["sd"] = par2, _["xi"] = par3);
        } else if (family == "Student-t") {
            Function dstd = fGarch["dstd"];
            m_dens = dstd(data_p, _["mean"] = par1, _["sd"] = par2, _["nu"] = par3);
        } else if (family == "Skew Student-t") {
            Function dsstd = fGarch["dsstd"];
            m_dens = dsstd(data_p, _["mean"] = par1, _["sd"] = par2, _["nu"] = par3, _["xi"] = par4);
        }
    }
    
    double obj = 0.0;
    for(int i = 0; i < n; i++) {
        double d = m_dens[i];
        if (std::isnan(d) || d <= 0) d = 1e-100;
        obj -= z_values[i] * std::log(d);
    }
    
    return obj;
}

// [[Rcpp::export]]
NumericVector eval_margin_cpp(NumericVector data_p, std::string family, NumericVector pars, std::string type) {
    double par1 = pars[0];
    double par2 = pars.size() > 1 ? pars[1] : 0.0;
    double par3 = pars.size() > 2 ? pars[2] : 0.0;
    double par4 = pars.size() > 3 ? pars[3] : 0.0;
    
    int n = data_p.size();
    NumericVector result(n);
    
    if (family == "Cauchy") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") result[i] = R::dcauchy(data_p[i], par1, par2, 0);
            else if (type == "cdf") result[i] = R::pcauchy(data_p[i], par1, par2, 1, 0);
            else if (type == "quant") result[i] = R::qcauchy(data_p[i], par1, par2, 1, 0);
        }
    } else if (family == "Gamma") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") result[i] = R::dgamma(data_p[i], par1, 1.0/par2, 0);
            else if (type == "cdf") result[i] = R::pgamma(data_p[i], par1, 1.0/par2, 1, 0);
            else if (type == "quant") result[i] = R::qgamma(data_p[i], par1, 1.0/par2, 1, 0);
        }
    } else if (family == "Logistic") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") result[i] = R::dlogis(data_p[i], par1, par2, 0);
            else if (type == "cdf") result[i] = R::plogis(data_p[i], par1, par2, 1, 0);
            else if (type == "quant") result[i] = R::qlogis(data_p[i], par1, par2, 1, 0);
        }
    } else if (family == "Loglogistic") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") {
                if (data_p[i] > 0) {
                    double term = std::pow(data_p[i] * par2, par1);
                    result[i] = (par1 * par2 * std::pow(data_p[i] * par2, par1 - 1.0)) / std::pow(1.0 + term, 2.0);
                } else result[i] = 0.0;
            } else if (type == "cdf") {
                if (data_p[i] > 0) {
                    double term = std::pow(data_p[i] * par2, par1);
                    result[i] = term / (1.0 + term);
                } else result[i] = 0.0;
            } else if (type == "quant") {
                result[i] = (1.0 / par2) * std::pow(data_p[i] / (1.0 - data_p[i]), 1.0 / par1);
            }
        }
    } else if (family == "Lognormal") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") result[i] = R::dlnorm(data_p[i], par1, par2, 0);
            else if (type == "cdf") result[i] = R::plnorm(data_p[i], par1, par2, 1, 0);
            else if (type == "quant") result[i] = R::qlnorm(data_p[i], par1, par2, 1, 0);
        }
    } else if (family == "Normal") {
        for(int i = 0; i < n; i++) {
            if (type == "pdf") result[i] = R::dnorm(data_p[i], par1, par2, 0);
            else if (type == "cdf") result[i] = R::pnorm(data_p[i], par1, par2, 1, 0);
            else if (type == "quant") result[i] = R::qnorm(data_p[i], par1, par2, 1, 0);
        }
    } else if (family == "Skew Normal" || family == "Student-t" || family == "Skew Student-t") {
        Environment fGarch = Environment::namespace_env("fGarch");
        if (family == "Skew Normal") {
            Function func = (type == "pdf") ? fGarch["dsnorm"] : ((type == "cdf") ? fGarch["psnorm"] : fGarch["qsnorm"]);
            result = func(data_p, _["mean"] = par1, _["sd"] = par2, _["xi"] = par3);
        } else if (family == "Student-t") {
            Function func = (type == "pdf") ? fGarch["dstd"] : ((type == "cdf") ? fGarch["pstd"] : fGarch["qstd"]);
            result = func(data_p, _["mean"] = par1, _["sd"] = par2, _["nu"] = par3);
        } else if (family == "Skew Student-t") {
            Function func = (type == "pdf") ? fGarch["dsstd"] : ((type == "cdf") ? fGarch["psstd"] : fGarch["qsstd"]);
            result = func(data_p, _["mean"] = par1, _["sd"] = par2, _["nu"] = par3, _["xi"] = par4);
        }
    }
    
    return result;
}

// [[Rcpp::export]]
NumericMatrix eval_all_margins_cpp(NumericMatrix data, CharacterVector families, NumericMatrix params, std::string type) {
    int n = data.nrow();
    int p = data.ncol();
    NumericMatrix result(n, p);
    
    for (int j = 0; j < p; j++) {
        NumericVector col_data = data( _, j);
        std::string family = as<std::string>(families[j]);
        NumericVector pars = params( _, j);
        result( _, j) = eval_margin_cpp(col_data, family, pars, type);
    }
    
    return result;
}
