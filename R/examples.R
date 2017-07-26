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

# Process development example, BH^2 (2nd ed) p. 200
pdEff = structure(
    c(-8, 24, 1, -2.25, 0.75, -1.25, -0.75, -5.5, 0, 4.5, 
      0.5, -0.25, -0.25, -0.75, -0.25), 
    .Names = c("C", "T", "CT", "P", "CP", "TP", "CTP", 
               "c", "Cc", "Tc", "CTc", "Pc", "CPc", "TPc", "CTPc"),
    mean = 72.25
)

# Simulated example
viseEff = structure(
    c(-6.88,  2.37,  0.25,  9.74,  3.76, -9.24, -5.59, 
      -0.49, -1.14,  9.66,  9.68, -7.59,  1.95, -2.71, -9.69),
    .Names = LETTERS[1:15],
    mean = 37.8
)

# BH^2 bike exper 2nd ed., p.245
bikeEff = structure(
    c(3.5, 12, 22.5, 1, .5, 1, 2.5),
    .Names = c("Seat", "Dyn", "Gear", "HBar", "Coat", "Bkfst", "Tires"),
    mean = 66.5
)

# Shrinkage effects BH^2 p.267 (w/ corrections, renamed, in std order of A, B, C, D)
shnkEff = structure(
    c(6.3375, 6.2125, 2.5875, -5.7375, 0.0375, 7.5125, 0.5375, 6.8875, 
      4.1625, 24.3875, 2.9125, 9.0875, 6.7125, -14.1625, 0.7125), 
    .Names = c("lTen", "lSpd", "mTmp", "lDie", "cMtl", "bTen", "Cool", 
               "lOD", "lTmp", "wTyp", "Scrn", "lMtl", "cDie", "wDia", "Spd"),
    mean = 29.10625
)

# effects for log(variance)
shnkDisp = structure(
    c(-0.8579, -0.8101, -0.2117, 0.4031, 0.1521, 1.1184, 0.672, 0.2932, 
      -0.2785, -0.1365, -0.9337, -0.2171, 0.0822, 0.7488, 0.6332), 
    .Names = c("lTen", "lSpd", "mTmp", "lDie", "cMtl", "bTen", "Cool", 
               "lOD", "lTmp", "wTyp", "Scrn", "lMtl", "cDie", "wDia", "Spd"),
    mean = 2.060015
)
