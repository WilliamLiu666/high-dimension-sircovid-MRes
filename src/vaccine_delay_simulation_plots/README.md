The structure of the task can be broken down into two primary functions found
in `plots.R`, with callback occurring in `script.R`. The code has been written
with a combination of ggplot2 and patchwork for speed and flexibility.

--------------------------------------------------------------------------------

IMPORTANT!!!

It is ESSENTIAL that you spend some time manually testing what values work for
the plot that you wish to plot. Due to the logic of the code, it is made to be
highly customisable as previous iterations that have attempted to automate plots
have not worked as expected. Therefore we have opted to trade speed for 
precision.

--------------------------------------------------------------------------------

The main function which generates these is the following;

ggplot_interventions <- function(State,
                                 row_facet,
                                 column_facet,
                                 rstats,
                                 school_closures,
                                 y_lim,
                                 force_scale,
                                 text_ylim = y_lim,
                                 x_offset,
                                 ylab_offset,
                                 plot_text,
                                 log)
                                 
State (string); the variable wished to be observed on the y-axes

row_facet (string); the parameter which will vary across columns

column_facet (string); the parameter which will vary across rows

rstats (output from function); this is the output of the function r_stats()
which generates daily r growth rates to plot as text (see `support.R`)

school_closuresl (boolean); if TRUE this will plot squre rectangles across
the x-axes that show the dates for which schools are closed, currently this is
manually defined wihin the function `ggplot_trajectory` but should be passed 
from a .csv (TODO). if FALSE it won't plot anything.

y_lim (vector of numeric); vector sets the maximum limit across the facet rows,
therefore you must ensure this is of length equal to the number of rows you wish
to plot.

force_scale (numeric); a numeric that acts as a force modifier to separate text
which is written on the plots. this modifier is taken in the function
`force_pull()` within `ggplot_interventions`.

text_ylim (vector of numeric); vector sets the "median" position of the text 
defined in the `r_stats()` output and must be the same length as y_lim.

x_offset (date); to be inputted as `.as.Date("20xx-xx-xx")`, this will set the
x-axis location of the text chunks that are plotted. As an example log-plots
will usually have an x_offset closer to the origin due to the decay trend. while
natural-plots are best placed at the end of the siumulation as the epidemic will
have ended and text won't overlap.

ylab_offset (numeric); this value shifts the location of the y-axes label

plot_text (boolean); if TRUE the plots will be returned with growth rates
plotted as text, if FALSE this will not be plotted.
recommended plot_text = FALSE if displaying Rt.

log (boolean); if TRUE the plot will transform the scale to log, if false
the scale will remain natural.

--------------------------------------------------------------------------------

If there are any issues with the plot then some manual input and debugging using
`browser()` (or your preferred method of debugging functions) will be needed to
be placed either in `ggplot_interventions` if dealing with aesthetics or if
dealing with the "core" plot then please explore `ggplot_trajectory` as this
deals directly with ggplot2 logic.