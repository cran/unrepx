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

### Various SE estimates ###

# Accessor to all _pse functions
PSE = function(effects, method = "Zahn", verbose = FALSE) {
    pse = get(paste0(method, "_pse"))
    if (!is.null(setup <- attr(pse, "setup"))) {
        parm = setup(length(effects))
        result = pse(effects, parm)
        if(verbose) {
            cat("Parameters used by PSE function:\n")
            print(parm)
        }
    }
    else 
        result = pse(effects)
    names(result) = paste0(method, "_PSE")
    result
}


# Scaled median - 1st step of Lenth PSE
SMedian_pse = function(effects)
    1.5 * stats::median(abs(effects))

# Root mean square
RMS_pse = function(effects) {
    sqrt(sum(effects^2))
}

# Lenth PSE
Lenth_pse = function(effects) {
    abseff = abs(effects)
    s0 = 1.5 * stats::median(abseff)
    1.5 * stats::median(abseff[abseff <= 2.5*s0])
}

# Dong93
Dong_pse = function(effects) {
    aeff = abs(effects)
    thresh = 3.75 * stats::median(aeff)
    sel = aeff[aeff <= thresh]
    sqrt(sum(sel^2) / length(sel))
}

# Daniel59
Daniel_pse = function(effects) {
    m = floor(.683*length(effects) + .5)
    sort(abs(effects))[m]
}

# Juan & Pena 1992
JuanPena_pse = function(effects) {
    abseff = sort(abs(effects))
    MAD = stats::median(abseff)
    top.n = length(effects)
    while (top.n != (top.n <- max(which(abseff <= 3.5 * MAD))))
        MAD = stats::median(abseff[seq_len(top.n)])
    MAD / .6578
}

# Zahn75 -- call Zahn.setup first
Zahn_pse = function(effects, parm) {
    abs.c = sort(abs(effects))[seq_len(parm$m)] 
    sum(parm$coef * abs.c)
}
attr(Zahn_pse, "setup") = function(n.effects) {
    m = floor(.683*n.effects + .5)
    q = (seq_len(m) - .375) / (n.effects + .25)
    z = stats::qnorm((1 + q)/2)
    list(m = m, coef = z / sum(z^2))
}

# Weighted version of Zahn
WZahn_pse = function(effects, parm) {
    abs.c = sort(abs(effects))[seq_len(parm$m)] 
    sum(parm$coef * abs.c)
}
attr(WZahn_pse, "setup") = function(n.effects) {
    m = floor(.683*n.effects + .5)
    q = (seq_len(m) - .375) / (n.effects + .25)
    z = stats::qnorm((1 + q)/2)
    w = m + .5 - seq_len(m)
    w[w > .65 * m] = .65 * m
    list(m = m, weight = w, coef = w*z / sum(w*z^2))
}



## Margin of error and SME
ME = function(effects, method = "Zahn", alpha = .05, ...) {
    pse = PSE(effects, method)
    rd = .getrefdist(length(effects), method, ...)
    c(ME = pse * quantile(rd$abst, 1 - alpha), 
      SME = pse * quantile(rd$max.abst, 1 - alpha))
}



# simulated reference dist of |t| and max |t|
# returns a list:
#    abst: dist of |t| under complete null (n x nsets matrix)
#    max.abst: dist on max(|t|) (vector of length nsets)
#    signature of form 'methodName_n.effects'
# If called with save = TRUE, result is saved in global variable .Last.ref.dist
#                                   update: in options("unrepx.lastrefdist")
ref.dist = function(method, n.effects, nsets, save = TRUE) {
    psefun = get(paste0(method, "_pse"))
    sig = paste(method, n.effects, sep="_")
    if(missing(nsets))
        nsets = ceiling(40000 / n.effects)
    X = matrix(rnorm(n.effects * nsets), nrow = n.effects) # each column is a sample
    # run setup code, if any
    if (!is.null(setup <- attr(psefun, "setup"))) {
        parm = setup(n.effects)
        abst = apply(X, 2, function(x) abs(x) / psefun(x, parm))
    }
    else
        abst = apply(X, 2, function(x) abs(x) / psefun(x))
    max.abst = apply(abst, 2, max)
    result = list(abst = abst, max.abst = max.abst, sig = sig)
    class(result) = "eff_refdist"
    if (save)
        options(unrepx.lastrefdist = result)
        #assign(".Last.ref.dist", result, envir = .GlobalEnv)
    result
}

print.eff_refdist = function(x, ...) {
    info = strsplit(x$sig, "_")[[1]]
    cat("Reference distribution of null effects\n")
    cat(paste0("Method: '", info[1], "', # effects = ", info[2], 
               ", # null samples = ", length(x$max.abst), "\n"))
    invisible()
}

# Get last ref.dist if it matches, else generate a new one
# May specify arguments in a list 'opts' rather than as arguments
.getrefdist = function(n.effects, method, nsets, save = TRUE, opts) {
    #refdist = try(get(".Last.ref.dist"), silent = TRUE)
    #if (inherits(refdist, "try-error")) 
    #    refdist = list(sig = "")
    refdist = getOption("unrepx.lastrefdist", list(sig = ""))
    sig = paste(method, n.effects, sep = "_")
    if(refdist$sig != sig) {
        if (!missing(opts)) {
            for (nm in names(opts))
                assign(nm, opts[[nm]])
        }
        refdist = ref.dist(method, n.effects, nsets = nsets, save = save)
    }
    refdist
}

eff.test = function(effects, method = "Zahn", pareto = TRUE, refdist, save = TRUE) {
    n.effects = length(effects)
    if(pareto) {
        ord = order(abs(effects))
        effects = effects[rev(ord)]
    }
    pse = PSE(effects, method)
    tval = effects / pse
    if(missing(refdist))
        refdist = .getrefdist(n.effects, method, save = save)
    pval = sapply(abs(tval), function(abst) mean(refdist$abst >= abst))
    spval = sapply(abs(tval), function(abst) mean(refdist$max.abst >= abst))
    result = data.frame(effect = effects, PSE = pse, t.ratio = round(effects/pse, 3),
               p.value = round(pval, 4), simult.pval = round(spval, 4))
    names(result)[2] = paste0(method, "_PSE")
    result
}


