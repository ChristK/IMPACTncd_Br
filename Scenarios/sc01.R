## Brazil Salt Policy model: a decision support tool for primary prevention of NCDs
## Copyright (C) 2019 Chris Kypridemos

## Brazil Sodium Policy model is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.


# This is the first policy scenario.
# The user needs to specify changes in added salt and salt from other sources
# all changes are additional to those in baseline scenario sc00 (net-effect)


## User inputs -----
sc_inputs <- list()
sc_inputs$names <- "sc01" # REMEMBER TO UPDATE for new scenarios

sc_inputs$year_change_starts      <- 2013 # Year the change in salt starts
sc_inputs$year_change_completed   <- 2017 # Year the change in salt is completed (linear diffusion). Change then remains for all consequent years
sc_inputs$men_salt_added_change   <- 0 # Mean change in added salt for men (g/d)
sc_inputs$women_salt_added_change <- 0 # Mean change in added salt for women (g/d)
sc_inputs$men_salt_other_change   <- -0.249  # Mean change in other salt for men (g/d)
sc_inputs$women_salt_other_change <- -0.249  # Mean change in other salt for women (g/d)

## Model logic. DO NOT ALTER ----
cat(paste0("Loading ", sc_inputs$names, "...\n"))
sc_inputs$names <- paste0("_", sc_inputs$names)
sc_inputs$year_change_starts <- sc_inputs$year_change_starts - 2000L
sc_inputs$year_change_completed <- sc_inputs$year_change_completed - 2000L

source(exprs = str2expression(gsub("____changeme",  sc_inputs$names,
                                   readLines(file.path("Scenarios", "template_policy_scenario.R")))), local = TRUE
)
