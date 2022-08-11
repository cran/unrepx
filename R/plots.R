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

### Plots for visual identification of effects


# Half-normal plot (or normal plot)
# if col = TRUE, we color pos blue, neg red, zero black. if FALSE, use black. Else use colors in col
hnplot = function(effects, ref = TRUE, half = TRUE, horiz = TRUE, method = "Zahn",
                  a = .375, col = half, pch = 16, ID = FALSE, alpha, ...) {
    n.effects = length(effects)
    abseff = abs(effects)
    q = (rank(abseff) - a) / (n.effects + 1 - 2*a)
    hns = qnorm((1 + q) / 2)
    
    if (is.logical(col)) {
        if (col)
            col = sapply(zapsmall(effects), 
                         function(x) ifelse(x < 0, "red", ifelse(x > 0, "blue", "black")))
        else
            col = "black"
    }

    if (half) {
        eff = abseff
        elab = "Absolute effects"
        score = hns
        slab = "Half-normal quantiles"
    }
    else {
        eff = effects
        elab = "Effects"
        score = qnorm((rank(eff) - a) / (length(eff) + 1 - 2*a))
        slab = "Normal quantiles"
    }
    
    if (ref)
        slope = PSE(effects, method)
    
    if (!missing(alpha)) {
        refdist = .getrefdist(n.effects, method, save = FALSE)
        MEs = slope * c(quantile(refdist$abst, 1 - alpha), quantile(refdist$max.abst, 1 - alpha))
        MElabs = c("ME", "SME")
        if (half)
            efflim = range(c(eff, MEs[1]))
        else {
            efflim = range(c(-MEs[1], eff, MEs[1]))
            MEs = c(-rev(MEs), MEs)
            MElabs = c("-SME", "-ME", "ME", "SME")
        }
    }
    else
        efflim = range(eff)
    
    if (is.numeric(ID)) {
        thresh = ID[1]
        ID = TRUE
    }
    else
        thresh = Inf
    
    if(ID) {
        if (is.null(names(abseff)))
            names(abseff) = seq_along(abseff)
        idx = which(abseff >= thresh)          # indices to auto-ID
        if (!is.infinite(thresh) && (length(idx) == 0))
            ID = FALSE
    }
        

    if(horiz) {
        plot(score ~ eff, xlim = efflim, ylab = slab, xlab = elab, col = col, pch = pch, ...)
        if(ref) 
            abline(0, 1 / slope)
        if (!missing(alpha)) {
            sapply(MEs, function(xx) abline(v = xx, lty = 2))
            axis(3, at = MEs, labels = MElabs, cex.axis = .75)
        }
        if (ID) {
            if (is.infinite(thresh)) 
                identify(eff, score, labels = names(abseff))
            else
                text(eff[idx], score[idx], names(abseff[idx]), 
                     pos = sapply(sign(eff[idx]), function(s) ifelse(s>0, 2, 4)))
        }
    }
    else {
        plot(eff ~ score, ylim = efflim, xlab = slab, ylab = elab, col = col, pch = pch, ...)
        if(ref) 
            abline(0, slope)
        if (!missing(alpha)) {
            sapply(MEs, function(xx) abline(h = xx, lty = 2))
            axis(4, at = MEs, labels = MElabs, cex.axis = .75)
        }
        if (ID) {
            if (is.infinite(thresh))
                identify(score, eff, labels = names(abseff))
            else
                text(score[idx], eff[idx], names(abseff[idx]),
                     pos = sapply(sign(eff[idx]), function(s) ifelse(s>0, 1, 3)))
        }
    }
    if (!half && ref) {
        abline(h = 0, lty = 3, col = "darkgray")
        abline(v = 0, lty = 3, col = "darkgray")
    }
    if (ref)
        title(sub = paste0("Reference SD = ", signif(slope, 4), " (", method, " method)"),
              cex.sub=.75, adj = 1)
    if(!missing(alpha)) {
        MEs = rev(MEs)
        title(sub = paste("ME =", signif(MEs[2], 3)," SME =", signif(MEs[1], 3)),
              cex.sub=.75, adj = 0)
    }
    
    invisible(list(eff = eff, abseff = abseff, hnscore = hns))
}



### Version with automatic bins
dot.plot = function(x, pch = 16, cex.dot = 1, spacing = 1, xlab, xlim = range(x), ...) {
    if(missing(xlab))
        xlab = as.character(substitute(x))
    
    # max(x) - min(x)]
    xspan = diff(range(x))
    
    # determine inverse ranks
    idx = seq_along(x)
    idx[order(x)] = idx
    
    # make a blank plot
    plot(x, rep(0, length(x)), type = "n", axes = FALSE, xlab = xlab, ylab = "", xlim = xlim)
    
    # draw scale
    axis(1)
    ylow = par("usr")[3]
    ytop = par("usr")[4]
    abline(h = ylow) # extend to full width
    
    # Just in case we ID them later
    if(is.null(names(x)))
        names(x) = seq_along(x)
    
    
    env = new.env()
    env$orig.x = x # used by dot.id()
    env$cex.dot = cex.dot
    
    RG = if(dev.interactive()) recordGraphics
    else function(expr, list, env) expr
    
    # draw points and support resizing in both directions
    RG({
        cxy = par("cxy") * cex.dot
        bins = as.integer(1.25 * xspan / cxy[1])
        xinc = diff(pretty(x, n = bins)[1:2])
        rx = xinc * round(x / xinc, 0)
        freq = table(rx)
        xx = rep(as.numeric(names(freq)), freq) [idx]
        yy = unlist(lapply(freq, seq_len)) [idx]
        yinc = env$yinc = 0.5 * spacing * cxy[2]
        yyv = ylow + yinc * (yy - .5)
        points(xx, yyv, pch = pch, cex = env$cex.dot, ...)
        clip = which(yyv > ytop)
        if(length(clip) > 0)
            text(xx[clip], rep(ytop + .5 * yinc, length(clip)), "(more)", col = "red", cex=.5, xpd = TRUE)
        # Additional labels user may add programatically later
        # by defining x.id, height, cex.id, col.id in env
        if(!is.null(env$x.id)) {
            rxl = xinc * round(env$x.id / xinc)
            for (xl in rxl) {
                labs = names(rx[rx == xl])
                yl = ylow + cxy[2] * (env$height.id - 1 + env$cex.id * seq_along(labs)) 
                text(xl, yl, labs, col = env$col.id, cex = env$cex.id)
            }
        }
    },
    list(), env)
    
    invisible(env)
}

# "identify" function for dot plots - specify environment returned by dot.plot
# with modify = TRUE, previously identified points are used but user can specify different height etc. 
dot.id = function(env, height.id = 2, cex.id = 1, col.id = "black") {
    if (!dev.interactive())
        stop("Graphics device is not interactive")
    ylow = par("usr")[3]
    x = env$orig.x
    pts = identify(x, rep(ylow, length(x)), pos = FALSE, plot = FALSE)
    if (length(pts) > 0) {
        modify = TRUE
        env$x.id = env$orig.x[pts]
        env$height.id = height.id
        env$cex.id = cex.id
        env$col.id = col.id
    }
    message("You may need to refresh or slightly resize the plot")
    invisible(env$x.id)
}

dot.mod = function(env, ...) {
    if (!dev.interactive())
        stop("Graphics device is not interactive")
    dots = list(...)
    names(dots) = sapply(names(dots), function(nm) match.arg(nm, c("cex.dot", "height.id", "col.id", "cex.id")))
    for (nm in names(dots))
        env[[nm]] = dots[[nm]]
    message("You may need to refresh or slightly resize the plot")
}

# Reference plot
refplot = function(effects, ref = TRUE, half = TRUE, method = "Zahn", 
                   col = half, guides = FALSE, ID = FALSE, pch = 16, xlab, xlim, ...) {
    
    n.effects = length(effects)
    if(missing(xlab))
        xlab = ifelse(half, "Absolute effects", "Effects")
    ae = abs(effects)

    if (is.logical(col)) {
        if (col)
            col = sapply(zapsmall(effects), 
                         function(x) ifelse(x < 0, "red", 
                                            ifelse(x > 0, "blue", "black")))
        else
            col = "black"
    }
    
    if (half)
        effects = abs(effects)
    xlim = range(c(0, effects))
    
    env = dot.plot(effects, xlim = xlim, xlab = xlab, pch = pch, col = col, ...)
    
    if(is.character(ref)) {
        curvetype = match.arg(ref, c("normal", "simulated"))
        ref = TRUE
        guides = FALSE
    }
    else
        curvetype = "normal"
        
    ytop = ifelse(guides, 1.23, 1)
    ylow = par("usr")[3]
    yscal = (par("usr")[4] - ylow) / ytop
    
    if(ref) {
        stderr = PSE(effects, method)
        x = seq(-3 * stderr, 3 * stderr, len = 81)
        if (half)
             x = x[41:81]
 
        if (curvetype == "normal") {
            lines(x, ylow + yscal * exp(-.5*(x / stderr)^2))
            lines(c(0, 0), ylow + yscal * c(0,1), lty = 2)
            title(sub = paste0("Reference curve SD = ", signif(stderr, 4), " (", method, " method)"),
                  cex.sub=.75, adj = 1)
        }
        else { # curvetype = "simulated"
            refdist = .getrefdist(n.effects, method)
            if(half)
                DE = stats::density(stderr * refdist$abst, bw = "SJ")
            else
                DE = stats::density(stderr * c(refdist$abst, -refdist$abst), bw = "SJ")
            lines(DE$x, ylow + yscal * (DE$y / max(DE$y)))
            title(sub = paste("Ref dist:", method, "method for", n.effects, "effects"),
                  cex.sub=.75, adj = 1)
        }
        if (guides) {
            e12 = exp(-.5)
            lines(c(-2*stderr, 0, 2*stderr), ylow + yscal * c(0, 2*e12, 0), lty = 3)
            points(c(-stderr, stderr), ylow + yscal * c(e12, e12), cex = .5)
            lines(c(0, 0), ylow + yscal * c(1, ytop), lty = 3)
        }
    }
    
    if (is.numeric(ID)) {
        thresh = ID[1]
        ID = TRUE
    }
    else
        thresh = Inf
    
    if (ID) {
        if (is.infinite(thresh))
            dot.id(env)
        else {
            env$x.id = effects[abs(effects) >= thresh]
            env$height.id = 2
            env$cex.id = 1
            env$col.id = "black"
        }
    }
    if(!dev.interactive()) 
        { # show the text since we won't be re-drawing
        text(env$x.id, ylow + yscal * env$height.id *env$yinc, names(env$x.id))
    }

    invisible(env)
}


# pareto plot
parplot = function (effects, pareto = TRUE, absolute = TRUE, horiz = FALSE, 
          col = absolute, critvals, method = "Zahn", 
          alpha = .05, refdist, sim.opts, ###save = TRUE, nsets, 
          ylab = "Estimated effects", 
          top = n.effects, cex.annot = .75, ...) 
{
    n.effects = length(effects)
    if (is.null(names(effects)))
        names(effects) = seq_along(effects)
    
    if (is.logical(col)) {
        if (col)
            col = sapply(effects, 
                         function(x) ifelse(x < 0, "pink", "lightblue"))
        else
            col = rep("lightgray", n.effects)
    }
    
    if (missing(critvals)) {
        if(missing(refdist))
            refdist = .getrefdist(n.effects, method, opts = sim.opts)
        critvals = c(quantile(refdist$abst, 1 - alpha),
                    quantile(refdist$max.abst, 1 - alpha))
        critmeth = method
        pse = PSE(effects, method)
    }
    else {
        critvals = c(critvals, NA)  # SME critval is NA if not provided
        critmeth = "Manual"
        pse = 1
    }
    cutoffs = pse * c(-critvals[2], -critvals[1], 0, 
                                 critvals[1], critvals[2])
    cutlabs = c("-SME", "-ME", "0", "ME", "SME")
    names(cutoffs) = cutlabs
    
    if (pareto || (top < n.effects)) {
        ord = order(abs(effects))
        effects = effects[rev(ord)][seq_len(top)]
        col = col[rev(ord)][seq_len(top)]
    }
    
    ylim = 1.05 * range(c(effects, cutoffs[c(2, 4)]))
    if (absolute) {
        effects = abs(effects)
        ylim[2] = max(abs(ylim))
        ylim[1] = 0
    }
    if (!horiz) {
        barplot(effects, ylim = ylim, ylab = ylab, col = col, 
                ...)
        abline(h = 0)
        axis(side = 4, at = cutoffs, labels = cutlabs)
        sapply(cutoffs[c(1, 2, 4, 5)], function(x) abline(h = x, 
                                                          lty = 2))
        axis(side = 2, at = 2 * ylim, labels = c("", ""))
    }
    else {
        barplot(rev(effects), horiz = TRUE, xlim = ylim, xlab = ylab, 
                col = rev(col), ...)
        abline(v = 0)
        axis(side = 3, at = cutoffs, labels = cutlabs)
        sapply(cutoffs[c(1, 2, 4, 5)], function(x) abline(v = x, 
                                                          lty = 2))
        axis(side = 1, at = 2 * ylim, labels = c("", ""))
    }
    melab = paste0(critmeth, " critical values: ME = ", signif(cutoffs[4], 3), 
                   ", SME = ", signif(cutoffs[5], 3))
    
    # function to place annotations
    annot = function(annot, adj, col = "black") {
        title(sub = annot, adj = adj, cex.sub = cex.annot, col.sub = col)
    }
    if(cex.annot > 0)
        annot(melab, 0)
    if(top < n.effects) {
        if(cex.annot <= 0) cex.annot = .75
        annot(paste("Only the largest", top, "effects shown"), 1, "blue")
    }
    invisible(cutoffs[4:5])
}


