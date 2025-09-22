suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(tidyverse)
    library(magrittr)
    library(glue)
})
###############
# Handling/formatting plots
make_ggtheme <- function(...){
    theme(
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            ),
        ...
    )
}

scale_y_axis <- function(
    figure,
    mode='',
    log_base=10,
    axis_label_accuracy=2,
    n_breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # Scale y axis based on mode argumetn
    if (mode == 'pct') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            # limits=limits,
            expand=expand,
            n.breaks=n_breaks,
            labels=label_percent(),
            ...
        )
    } else if (mode == 'mb') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_continuous(
            # limits=limits,
            expand=expand,
            n.breaks=n_breaks,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis_label_accuracy
                ),
            ...
        )
    } else if (mode == 'log10') {
        figure +
        coord_cartesian(ylim=limits) +
        scale_y_log10(
            # limits=limits,
            expand=expand,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log_base,
                    digits=axis_label_accuracy,
                    signed=FALSE
                ),
            ...
        )
    } else if (mode == '') {
        if (is.null(limits)) {
            figure + 
            scale_y_continuous(
                expand=expand,
                labels=
                    function(x) {
                        format(x, digits=axis_label_accuracy)
                    }
            )
        } else {
            figure + 
            coord_cartesian(ylim=limits) +
            scale_y_continuous(
                # limits=limits,
                expand=expand,
                ...
            )
        }
        # coord_cartesian(ylim=limits) +
        # scale_y_continuous(...)
    }
}

add_faceting <- function(
    figure,
    facet_group=NULL,
    facet_col=NULL,
    facet_row=NULL,
    ...){
    # Facet as specified
    if (!is.null(facet_col) & !is.null(facet_row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                cols=vars(!!sym(facet_col)),
                ...
            )
    } else if (!is.null(facet_row)) {
        figure <- 
            figure +
            facet_grid2(
                rows=vars(!!sym(facet_row)),
                ...
            )
    } else if (!is.null(facet_col)) {
        figure <- 
            figure +
            facet_grid2(
                cols=vars(!!sym(facet_col)),
                ...
            )
    } else if (!is.null(facet_group)) {
        figure <- 
            figure +
            facet_wrap2(
                vars(!!sym(facet_group)),
                ...
            )
    }
    figure
}

post_process_plot <- function(
    figure,
    theme_obj=NULL,
    # add_theme=TRUE,
    facet_row=NULL,
    facet_col=NULL,
    facet_group=NULL,
    scales='fixed',
    scale_mode='',
    log_base=10,
    axis_label_accuracy=2,
    n_breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    # legend.position='right',
    # legend.ncols=1,
    # axis.text.x=element_text(angle=35, hjust=1),
    ...){
    figure <- 
        add_faceting(
            figure,
            facet_row=facet_row,
            facet_col=facet_col,
            facet_group=facet_group,
            scales=scales
        )
    # Set y-axis scaling (log, Mb, percent etc.)
    figure <- 
        figure %>% 
        scale_y_axis(
            mode=scale_mode,
            log_base=log_base,
            axis_label_accuracy=axis_label_accuracy,
            n_breaks=n_breaks,
            limits=limits,
            expand=expand
        )
    # Add theme elements, as either an object or individual args
    if (!is.null(theme_obj)) {
        figure <- figure + theme_obj
    } else {
        figure <- figure + make_ggtheme(...)
    } 
    figure
}

###############
# Make tabs per plot in Rmd
plot_figure_tabs <- function(
    plot.df,
    group_col,
    plot_fnc,
    header_lvl,
    nl_delim,
    return_figure=FALSE,
    ...){
    # message(paste(header_lvl, group_col, collapse=','))
    # print(table(plot.df[[group_col]]))
    plot.df[[group_col]] %>% 
    as.factor() %>% 
    droplevels() %>% 
    levels() %>%
    sapply(
        function(group_value, plot.df, plot_fnc, header_lvl, group_col, nl_delim, return_figure){
            figure <- 
                plot.df %>%
                # filter(get({{group_col}}) == group_value) %>%
                filter(!!sym(group_col) == group_value) %>%
                plot_fnc(...)
                if (return_figure) {
                    figure
                } else {
                cat(
                    strrep('#', header_lvl), group_value,
                    # nl_delim, "Rows per df", nrow(plot.df),
                    nl_delim
                )
                print(figure)
                cat(nl_delim)
            }
        },
        plot.df=plot.df,
        plot_fnc=plot_fnc,
        header_lvl=header_lvl,
        group_col=group_col,
        nl_delim=nl_delim,
        return_figure=return_figure,
        simplify=FALSE,
        USE.NAMES=TRUE
    )
}

make_tabs_recursive <- function(
    plot.df, 
    group_cols,
    current_header_lvl,
    plot_fnc,
    tabset_format,
    nl_delim,
    return_figure=FALSE,
    ...){
    # cat("LENGTH OF GROUP COLS", length(group_cols), group_cols, "\n\n\n")
    if (length(group_cols) == 1) {
        plot_figure_tabs(
            plot.df=plot.df, 
            group_col=group_cols[1],
            header_lvl=current_header_lvl,
            plot_fnc=plot_fnc,
            nl_delim=nl_delim,
            return_figure=return_figure,
            ...
        )
    } else {
        group_col <- group_cols[1]
        # message(paste(current_header_lvl, group_col, collapse=','))
        group_values <- 
            plot.df[[group_col]] %>% 
            as.factor() %>% 
            droplevels() %>% 
            levels()
        for (group_value in group_values) {
            cat(
                strrep('#', current_header_lvl), group_value, tabset_format,
                # nl_delim, "Rows:", nrow(plot.df), 
                nl_delim
            )
            make_tabs_recursive(
                plot.df=plot.df %>% filter(get({{group_col}}) == group_value),
                group_cols=group_cols[2:length(group_cols)],
                current_header_lvl=current_header_lvl + 1,
                plot_fnc=plot_fnc,
                tabset_format=tabset_format,
                nl_delim=nl_delim,
                ...
            )
        }
    }
}

make_nested_plot_tabs <- function(
    plot.df,
    group_cols,
    plot_fnc,
    max_header_lvl=2,
    tabset_format="{.tabset .tabset-pills}",
    nl_delim="\n\n\n",
    return_figure=FALSE,
    ...){
    cat(nl_delim)
    cat(strrep('#', max_header_lvl), tabset_format, nl_delim)
    plot.df %>% 
    make_tabs_recursive(
        group_cols=group_cols,
        current_header_lvl=max_header_lvl+1,
        plot_fnc=plot_fnc,
        tabset_format=tabset_format,
        nl_delim=nl_delim,
        return_figure=return_figure,
        ...
    )
    cat(nl_delim)
}

###############
# Basic Plots
plot_distribution <- function(
    plot.df,
    x_var='',
    fill_var=NULL, 
    label.size=3,
    label.y=2000,
    alpha=0.8,
    linewidth=0.75,
    x_off=0.075,
    binwidth=0.02,
    txt.angle=90,
    ...){
    # facet_group=NULL,
    # facet_col=NULL,
    # facet_row=NULL,
    # scales='fixed',
    # scale_mode='',
    # legend.position='right',
    # axis.text.x=element_text(angle=35, hjust=1),
    # Set fill group if specified
        # x_var='frequency'; facet_row='direction'; facet_col='significance'; alpha=0.8
    figure <- 
        if (is.null(fill_var)) {
            ggplot(
                plot.df,
                aes(x=.data[[x_var]])
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x_var]],
                    fill=.data[[fill_var]]
                )
            )
        }
    # make it a boxplot 
    figure <- 
        figure + 
        geom_histogram(
            aes(y=after_stat(count)),
            binwidth=binwidth,
            position='identity',
            alpha=alpha
        ) + 
        labs(y='# of Genes')
    # Handle faceting + scaling + theme options
    figure <- 
        figure %>% 
        post_process_plot(
            # facet_row=facet_row,
            # facet_col=facet_col,
            # facet_group=facet_group,
            # scales=scales,
            # scale_mode=scale_mode,
            # axis_label_accuracy=axis_label_accuracy,
            # legend.position=legend.position,
            # axis.text.x=axis.text.x,
            ...
        )
    taf1.df <- 
        plot.df %>%
        filter(EnsemblID == TAF1_ENSEMBLID) %>%
        mutate(TAF1.label=paste('TAF1', round(!!sym(x_var), 3), sep=': ')) %>% 
        mutate(
            x_offset=
                case_when(
                    !!sym(x_var) > max(!!sym(x_var)) / 2 ~ !!sym(x_var) - x_off,
                    TRUE ~ !!sym(x_var) - x_off
                )
        )
        # mutate(TAF1.label=glue('TAF1: {round(!!sym(x_var), 3)}'))
        # mutate(TAF1.label=glue('TAF1: {round({x_var}, 3)}'))
    taf1.figure <- 
        if (is.null(fill_var)) {
            figure +
            geom_vline(
                data=taf1.df,
                aes(xintercept=.data[[x_var]]),
                linewidth=linewidth,
                linetype='dashed'
            ) +
            geom_text(
                data=taf1.df,
                aes(
                    x=x_offset,
                    y=label.y,
                    label=TAF1.label
                ),
                angle=txt.angle,
                hjust=-0.25,
                size=label.size
             )
        } else {
            figure +
            geom_vline(
                data=taf1.df,
                aes(
                    xintercept=.data[[x_var]],
                    color=.data[[fill_var]]
                ),
                linewidth=linewidth,
                linetype='dashed'
            ) +
            geom_text(
                data=taf1.df,
                aes(
                    x=x_offset,
                    y=label.y,
                    label=TAF1.label,
                    color=.data[[fill_var]]
                ),
                angle=txt.angle,
                hjust=-0.25,
                size=label.size
             )
        }
    taf1.figure
}

plot_compare_results <- function(
    plot.df,
    color_var,
    full.BH.thresh=0.01,
    size=0.15,
    alpha=0.6,
    linewidth=0.01){
    figure <- 
        plot.df %>% 
        ggplot(
            aes(
                x=frequency,
                y=-log10(full.model.pvalue)
                # aes(color=.data[[color_var]])
                # color=full.model.significance
            )
        ) +
        geom_point(
            aes(color=.data[[color_var]]),
            size=size,
            alpha=alpha
        ) +
        geom_smooth(method=lm) +
        geom_hline(
            yintercept=-log10(full.BH.thresh),
            linetype='dashed',
            linewidth=linewidth
        ) +
        scale_y_continuous(
            limits=c(0, NA),
            expand=c(0.00, 0.00, 0.00, 0.00)
        ) +
        facet_grid(
            cols=vars(direction),
            rows=vars(significance),
            scales='fixed'
        ) + 
        theme(legend.position='top')
    # Color stuff
    if (color_var == 'full.model.direction') {
        figure <- 
            figure +
            scale_color_manual(
                values=
                    c(
                        "N.S"='#d3d3d3',
                        'Up in Full Model'='#ffa07a',
                        'Down in Full Model'='#add8e6'
                    ) 
            )
    } else if (color_var == 'log2FoldChange') {
        figure <- 
            figure +
            scale_color_gradient2(
                low='#add8e6',
                mid='#d3d3d3',
                high='#ffa07a'
            )
    }
    figure
}

