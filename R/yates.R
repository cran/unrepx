##############################################################################
#    Copyright (c) 2017 Russell V. Lenth                                     #
#                                                                            #
#    This file is part of the unrepx package for R (*unrepx*)                #
#                                                                            #
#    *unrepx* is free software: you can redistribute it and/or modify        #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *unrepx* is distributed in the hope that it will be useful,             #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with R and *unrepx*.  If not, see                                 #
#    <https://www.r-project.org/Licenses/> and/or                            #
#    <http://www.gnu.org/licenses/>.                                         #
##############################################################################

# Yates's algorithm

# Yates's (or reverse Yates's) algorithm for 2-level factors, y in standard order
#    y is vector of length 2^k; returns 'contrasts' (doubled regr coefs)
#    If y is vector of length 2^k - 1, optionally with "mean" attribute, returns predictions

yates =  function(y, labels = LETTERS, sep = "") {
    n = length(y)
    k = round(log(n) / log(2))
    if (round((nn <- 2^k) - n) != 0) {
        if (nn != n + 1)
            stop ("length(y) must be a power of 2")
        # else y is a set of effects...
        n = nn
        labels = rep("+", 10)
        lowlab = "-"
        reverse = TRUE
        if(is.null(ybar <- attr(y, "mean")))
            ybar = 0
        y = c(ybar, y/2)
        coef = c(-1, 1, 1, 1)
    }
    else {
        reverse = FALSE
        coef = c(1, 1, 1, -1) / 2
        lowlab = ""
    }
    names(y) = rep("", n)
    even = seq(2, n, by = 2)
    odd = even - 1
    for (i in 1:k) {
        y = c(coef[1]*y[even] + coef[2]*y[odd], coef[3]*y[even] + coef[4]*y[odd])
        names(y) = paste(names(y), rep(c(lowlab, labels[i]), each = n/2), sep = sep)
    }
    if (!reverse) {
        ybar = y[1]
        y = 2 * y[-1]
        attr(y, "mean") = ybar
    }
    y
}

# Generalized [reverse] Yates's algorithm. 
# Required argument nlevels, product of which must equal length(y) or length(y) - 1
# If the latter is true, the reverse Yates's algorithm is applied
gyates = function(y, nlevels, basis = "poly") {
    labs = as.character(seq_len(10))
    if (!is.null(nlev <- attr(y, "nlevels")))
        nlevels = nlev
    if (!is.null(bas <- attr(y, "basis")))
        basis = bas
    if ((n <- prod(nlevels)) != length(y)) {
        if (n != length(y) + 1)
            stop("Product of 'nlevels' must equal length of 'y'")
        else {      # length(y) == n - 1
            reverse = TRUE
            if (is.null(ybar <- attr(y, "mean")))
                ybar = 0
            y = c(0, y)
        }
    }
    else {      # length(y) == n
        reverse = FALSE
        ybar = mean(y)
        labs = c(".", labs)
    }
    
    nm = rep("", n)
    for (i in seq_along(nlevels)) {
        k = nlevels[i]
        m = matrix(y, nrow = k)
        bas = paste0(basis[1 + (i - 1) %% length(basis)], "_gyb")
        X = get(bas)(k)
        if (!reverse)
            X = t(X)
        y = as.numeric(t(X %*% m))
        if(i > 1)
            nm = as.character(t(matrix(nm, nrow = k)))
        nm = paste0(nm, rep(labs[1:k], each = n / k))
    }
    names(y) = nm
    if (reverse)
        y = y - mean(y) + ybar
    else {
        y = y[-1]
        attr(y, "mean") = ybar
        attr(y, "nlevels") = nlevels
        attr(y, "basis") = basis
    }
    y
}

# Provided orthonormal basis functions - user may add others
# 1st column should be 1/sqrt(k), others orthonormal contrasts
poly_gyb = function(k) {
    cbind(1/sqrt(k), contr.poly(seq_len(k)))
}

helmert_gyb = function(k) {
    X = cbind(1, contr.helmert(seq_len(k)))
    apply(X, 2, function(x) x / sqrt(sum(x^2)))
}

