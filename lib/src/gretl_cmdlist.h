/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* commands recognized by gretl */

#ifndef CMDLIST_H
#define CMDLIST_H

char *gretl_commands[] = {
    ";",     
    "add",
    "addto",
    "adf",      
    "ar",       
    "arch",
    "boxplot",
    "chow",     
    "coeffsum",
    "coint",
    "coint2",
    "corc", 
    "corr",     
    "corrgm",   
    "criteria",
    "critical",
    "cusum",
    "data",
    "delete",
    "diff",
    "else",
    "end",
    "endif",
    "endloop",
    "eqnprint",
    "equation",
    "fcast", 
    "fcasterr",
    "fit", 
    "freq",     
    "genr",     
    "gnuplot",  
    "graph", 
    "hausman",
    "hccm", 
    "help",    
    "hilu",    
    "hsk", 
    "if",
    "import",
    "info", 
    "label",
    "labels",
    "lad",
    "lags",    
    "ldiff",
    "leverage",
    "lmtest", 
    "logit",
    "logs",
    "loop",
    "meantest",
    "mpols",
    "multiply",
    "mvavg",
    "nls",
    "noecho",
    "nulldata",
    "ols",     
    "omit",
    "omitfrom",
    "open",
    "outfile",
    "panel",
    "pca",
    "pergm",
    "plot",    
    "pooled",
    "print",  
    "probit",
    "pvalue",  
    "quit", 
    "rename",
    "reset",
    "rhodiff", 
    "rmplot",
    "run",
    "runs",
    "scatters",
    "seed",
    "setobs",
    "setmiss",
    "shell",   
    "sim",     
    "smpl",
    "spearman",
    "square",  
    "store",   
    "summary",
    "system",
    "tabprint",
    "testuhat",
    "tsls",    
    "var",
    "varlist",
    "vartest",
    "wls"
}; 

#endif /* CMDLIST_H */
