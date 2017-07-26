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

### Routine to make "Daniel paper"

# n = number of effects to plot
# lin.div = number of divisions on linear scale (1 + lin.div lines plotted)
# ranks are centered using (i - a) / (n + 1 - 2*a)
# ranks says which axis has the ranks
daniel.paper = function(n.effects = 15, linear = 0:40, a = .375, 
                        rank.axis = c("y", "x"), 
                        type = c("half.normal", "normal"), 
                        theme = "blue", lin.lab, rank.lab) 
{
    theme = get(paste0(theme, "_thm"))
    rank.axis = match.arg(rank.axis)
    type = match.arg(type)
    x = linear
    
    ranks = seq_len(n.effects)
    qtl = (ranks - a) / (n.effects + 1 - 2*a)
    y = scores = switch(type, 
               half.normal = stats::qnorm((1 + qtl)/2),
               normal = stats::qnorm(qtl))
    
    major = theme$major
    minor = theme$minor
    axis = if(is.null(theme$axis)) major else theme$axis
        
    mod.major = c(which(linear %% 5 == 0)[1], 0) # index for starting major ticks
    
    if(missing(lin.lab)) lin.lab = ifelse(type=="normal", "Effects", "Absolute effects")
    if(missing(rank.lab)) rank.lab = "Rank"
    xlab = lin.lab
    ylab = rank.lab
    if (rank.axis == "x") { # swap axes and labels
        tmp = x; x = y; y = tmp
        tmp = xlab; xlab = ylab; ylab = tmp
        mod.major = rev(mod.major)
    }
    plot(range(x), range(y), type = "n", 
         xlab = xlab, ylab = ylab, axes = FALSE, col.lab = axis$col)
    par(usr = c(range(x), range(y))) # makes axis ends same as plot-region ends
    for (i in seq_along(x)) {
        if (i %% 5 == mod.major[1]) abline(v = x[i], col = major$col, lty = major$lty, lwd = major$lwd)
        else abline(v = x[i], col = minor$col, lty = minor$lty)
    }
    for (i in seq_along(y)) {
        if (i %% 5 == mod.major[2]) abline(h = y[i], col = major$col, lty = major$lty, lwd = major$lwd)
        else abline(h = y[i], col = minor$col, lty = minor$lty)
    }
    side = switch(rank.axis, y = 2, x = 1)
    sel = which(ranks %% 5 == 0)
    axis(side, at = scores[sel], labels = ranks[sel], 
         lwd = 0, lwd.ticks = axis$lwd, col.ticks = axis$col, col.axis = axis$col)

    invisible()
}

blue_thm = list(major = list(col = "blue", lty = 1, lwd = 1), 
                minor = list(col = "cyan", lty = 1))

bw_thm = list(major = list(col = "black", lty = 1, lwd = 1),
              minor = list(col = "black", lty = 3))

gray_thm = list(major = list(col = "black", lty = 1, lwd = 1),
                minor = list(col = "gray", lty = 1))

bold_thm = list(major = list(col = "black", lty = 1, lwd = 2),
                minor = list(col = "black", lty = 1, lwd = 1))

