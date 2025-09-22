###############
# Dependencies
TAF1_ENSEMBLID <- 'ENSG00000147133'
suppressPackageStartupMessages({
    library(sva)
    library(DESeq2)
    library(tictoc)
    library(magrittr)
    library(glue)
    library(purrr)
    library(tidyverse)
    library(devtools) # for session_info()
})
###############
# Misc
install_dependencies <- function(){
    install.packages(
        c(
            'devtools',
            'furrr',
            'future',
            'glue',
            'here',
            'magrittr',
            'purrr',
            'tictoc',
            'tidyverse'
        )
    )
    install.packages('biconductor')
    BiocManager::install(
        c(
            'DESeq2',
            'sva'
        )
    )
}

check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    results_fnc,
    ...){
    # Set read/write functions based on filetype
    output_filetype <- results_file %>% str_extract('\\.[^\\.]*$')
    if (output_filetype == '.rds') {
        load_fnc <- readRDS
        save_fnc <- saveRDS
    } else if (output_filetype %in% c('.txt', '.tsv')) {
        load_fnc <- read_tsv
        save_fnc <- write_tsv
    } else {
        stop(glue('Invalid file extesion: {extension}'))
    }
    # Now check if the results file exists and load it
    tic()
    if (file.exists(results_file) & !force_redo) {
        message(glue('{results_file} exists, not recomputing results'))
        if (return_data) {
            message('Loading cached results...')
            results <- load_fnc(results_file)
        }
    } else {
        if (file.exists(results_file) & force_redo) {
            message(glue('{results_file} exists, recomputing results anyways'))
        } else {
            message(glue('No cached results, generating: {results_file}'))
        }
        results <- results_fnc(...) %T>% save_fnc(results_file)
        # If results dont exist or force_redo is TRUE compute + cache results
        # Assumes save_fnc is of fomr save_fnc(result_object, filename)
    }
    toc()
    if (return_data) return(results) else return(invisible(NULL))
}

make_matrix_tidy <- function(
    counts.matrix,
    sample.metadata,
    row_names,
    join_col,
    value_names){
    # Make data tidy i.e. row is gene expression for 1 gene in 1 sample + metadata
    counts.matrix %>%
    as.data.frame() %>%
    rownames_to_column(var=row_names) %>%
    as_tibble() %>% 
    tidyr::pivot_longer(
        # where(is.numeric),
        cols=-all_of(c(row_names)),
        names_to=join_col,
        values_to=value_names,
    ) %>%
    # Add metadata for grouping samples
    left_join(
        sample.metadata,
        by=join_col
    )
}
###############
# Load model data
load_sample_metadata <- function(){
    SAMPLE_METADATA_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    select(
        SampleID,
        Region,
        Genotype,
        Brain,
        RIN,
        Library_Prep,
        Seq_Lane,
        Gender,
        Uniqmappedreads,
        Uniqmappedperc,
        Unmapped_mismatchperc,
        Unmapped_shortperc,
        `Exonic Rate`,
        `Intronic Rate`,
        `Intergenic Rate`,
        `Genes Detected`,
        `Mean 3' bias`,
        Exclude
    ) %>% 
    filter(!Exclude)
    # select(
    #     SampleID,
    #     # Brain,
    #     # CCXDP,
    #     Region,
    #     Genotype,
    #     # Genotype_Code,
    #     Gender,
    #     RIN,
    #     Library_Prep,
    #     Seq_Lane,
    #     Year_Collected,
    #     Batch_Brains,
    #     Uniqmappedreads,
    #     Uniqmappedperc,
    #     Unmapped_mismatchperc,
    #     Unmapped_shortperc,
    #     Exonic_Rate,
    #     Intronic_Rate,
    #     Intergenic_Rate,
    #     Genes_Detected,
    #     Mean_3p_bias,
    #     Median_ExonCV,
    #     TAF1_tpm,
    #     RptSize,
    #     MSH3,
    #     AAO,
    #     AAD,
    #     DisDur,
    #     PMI
    # ) %>%
    # {.}
    # rename(
    #     'reads.mapped.uniq'=Uniqmappedreads,
    #     'pct.mapped.uniq'=Uniqmappedperc,
    #     'pct.unmapped.mismatch'=Unmapped_mismatchperc,
    #     'pct.unmapped.short'=Unmapped_shortperc,
    #     'rate.exonic'=Exonic_Rate,
    #     'rate.intronic'=Intronic_Rate,
    #     'rate.intergenic'=Intergenic_Rate,
    #     'mean.3prime.bias'=Mean_3p_bias,
    #     'median.exonic.CV'=Median_ExonCV,
    #     'TPM.TAF1'=TAF1_tpm,
    #     # 'RepeatSize'=RptSize,
    #     # MSH3,
    #     # 'AgeAtOnset'=AAO,
    #     # 'AgeAtDeath'=AAD,
    #     # DisDur,
    #     # PMI
    # ) %>% 
    # rename_with(~ gsub('_', '.', .x)) %>%
    # mutate(
    #     DisDur=
    #         case_when(
    #             is.na(DisDur) ~ -Inf,
    #             TRUE ~ DisDur
    #         ),
    #     across(
    #         c(
    #             RptSize,
    #             AAO 
    #         ),
    #         ~ case_when(
    #             is.na(.x) ~ 0,
    #             TRUE ~ .x
    #         )
    #     ),
    #     log.RptSize=log(RptSize, base=2),
    #     across(
    #         c(
    #             Region,
    #             Genotype,
    #             Gender,
    #             Library_Prep,
    #             Seq_Lane,
    #             Year_Collected,
    #             Batch_Brains
    #         ),
    #         as.factor
    #     )
    # )
}


load_counts <- function(){
    RAW_COUNTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    column_to_rownames(var='EnsemblID') %>% 
    as.data.frame()
}

validate_models <- function(
    intercept,
    filter.variable,
    effect.variable,
    covariate.variables,
    sample.metadata,
    ...){
    # raise error if model arguments are invalid for building a DESeq2 model
    if (!(intercept %in% c(0, 1))) {
        stop('Invalid intercept, must be 0 or 1')
    }
    columns_undetected <- 
        c(
            filter.variable,
            effect.variable,
            covariate.variables
        ) %>% 
        unlist() %>% 
        unique() %>% 
        setdiff(colnames(sample.metadata)) %>%
        {.[!is.na(.) & !is.null(.)]}
    if (length(columns_undetected) != 0) {
        message(columns_undetected)
        stop('Columns argument specified missing in metadata, double check')
    }
}

load_model_specifications <- function(
    filepaths=NULL,
    output_dir=BOOTSTRAP_DIR,
    sample.metadata){
    {
        if (is.null(filepaths)) {
            MODEL_SPECIFICATIONS_DIR %>% 
            list.files(
                pattern='*.model',
                recursive=FALSE,
                full.names=TRUE
            )
        } else {
            filepaths
        }
    } %>% 
    read_delim(
        id='model.filename',
        delim='=',
        show_col_types=FALSE,
        col_names=FALSE
    ) %>%
    pivot_wider(
        id_cols=model.filename,
        names_from=X1,
        values_from=X2
    ) %>% 
    # {supp(readr::type_convert(.))} %>% 
    readr::type_convert(.) %>% 
    rowwise() %>% 
    mutate(covariate.variables=str_split(covariate.variables, pattern=',')) %T>% 
    pmap(
        .l=.,
        .f=validate_models,
        sample.metadata=sample.metadata
    ) %>% 
    ungroup() %>% 
    mutate(
        model.dir=
            file.path(
                output_dir,
                model.name
            ),
        info.filepath=
            file.path(
                model.dir,
                'bootstrap.parameters.txt'
            ),
        list.dir=
            file.path(
                model.dir,
                'bootstrap.samples.lists'
            ),
        results.dir=
            file.path(
                model.dir,
                'DEG.results'
            )
    )
}

load_bootstrap_chunks <- function(chunk.files){
    chunk.files %>% 
    tibble(chunk.file=.) %>% 
    rowwise() %>% 
    mutate(
        chunk_id=str_remove(basename(chunk.file), '.tsv'),
        bootstrap.chunk=
            read_tsv(
                chunk.file,
                show_col_types=FALSE
            ) %>%
            mutate(
                across(
                    ends_with('_samples'),
                    ~ str_split(.x, pattern=',')
                )
            ) %>%
            list()
    )
}
###############
# DESeq2 modelling
filter_genes_by_CPM <- function(
    sample_metadata,
    counts_matrix,
    sampleID.variable,
    gene_id_column,
    lib_sizes=NULL,
    CPM.cutoff=0.1,
    min.frac.samples.expressed=0.5) {
    # Pick genes to keep that are meaningfully expressed in >= 1 sample group
    # counts.matrix: raw expression matrix (genes x samples)
    # lib_sizes: library size for each sample
    # sample_conditions: a tibble with any number of categorical variables
    #                    matching to columns of counts.matrix
    # CPM.cutoff: cpm threshold for each expression value 
    # min.frac.samples.expressed: fraction of samples that have counts above CPM.cutoff
    # if (is.null(lib_sizes)) lib_sizes <- rowSums(counts.matrix)
    if (is.null(lib_sizes)) lib_sizes <- colSums(counts_matrix)
    # Calculate CPMs
    counts_matrix %>%
    {1000000 * (. / lib_sizes)} %>%
    # Pivot longer for summarizing
    make_matrix_tidy(
        sample.metadata=sample_metadata, 
        row_names=gene_id_column,
        join_col=sampleID.variable,
        value_names='counts'
    ) %>%
    # Group samples for cpm thresholding
    group_by(
        across(
            all_of(
                grep(
                    glue('{sampleID.variable}|counts'),
                    colnames(.),
                    value=TRUE,
                    invert=TRUE
                )
             )
        )
    ) %>%
    # Calculate fraction of samples > cutoff for each gene per group
    summarize(frac_above_cutoff=sum(counts > CPM.cutoff) / n()) %>% 
    # Keep any gene with enough expressed samples in at least 1 group
    summarize(to_keep=sum(frac_above_cutoff > min.frac.samples.expressed) > 0) %>%
    filter(to_keep) %>%
    pull(!!sym(gene_id_column)) 
}

check_covariates_valid <- function(
    sample_metadata,
    effect.variable,
    covariate.variables) {
    # Check if any variables are uniform (and thus pointeless to include)
    if (is.null(covariate.variables)) {
        return(list())
    } else if (all(is.na(covariate.variables))) {
        return(list())
    } else if (length(covariate.variables) == 0) {
        return(list())
    }
    # non-uniform variables
    non.uniform.variables <- 
        sapply(
            covariate.variables,
            function(covariate) length(unique(sample_metadata[[covariate]]))
        ) %>%
        {.[. > 1]} %>%
        names()
    # Full model matrix
    full.mod.matrix <- 
        c(effect.variable, non.uniform.variables) %>%
        paste(collapse="+") %>% 
        sprintf('~%s', .) %>% 
        formula() %>% 
        model.matrix(sample_metadata)
    # Now check which if any covariates are co-linear (can only include 1)
    # Matrix rank for each column if it was removed from here 
    # https://stats.stackexchange.com/questions/16327/testing-for-linear-dependence-among-the-columns-of-a-matrix
    colinearities <- 
        full.mod.matrix %>% 
        # get rank of matrix after removing each model column (covariate level)
        # All variables that share the max rank are colinear i.e. redundant 
        {
            sapply(
                colnames(.),
                function(colname) {
                    .[,!(colnames(.) == colname)] %>% 
                    qr() %>%
                    {.$rank}
                }
            )
        } %>%
        enframe(name='model.column', value='qr.rank') %>% 
        # Get variable names from model columns and how many groups they have
        rowwise() %>% 
        mutate(
            covariate= 
                str_extract(
                    model.column,
                    paste(c(effect.variable, covariate.variables), collapse='|')
                ),
        ) %>%
        mutate(
            is.effect=(covariate == effect.variable),
            levels=
                ifelse(
                    is.na(covariate),
                    0,
                    length(levels(sample_metadata[[covariate]]))
            )
        ) %>%
        ungroup() %>% 
        group_by(covariate) %>%
        mutate(qr.rank=max(qr.rank)) %>% 
        ungroup()
    if (length(unique(colinearities$qr.rank)) == 1) {
        return(list())
    }
    # full.mod.matrix
    # pick all variables not colinear with anything
    non.colinear.variables <- 
        colinearities %>% 
        filter(!is.effect, qr.rank < max(qr.rank)) %>%
        pull(covariate) %>% 
        unique()
    # pick one of the colinear variables with the most levels
    colinear.variables <- 
        colinearities %>% 
        filter(!is.na(covariate)) %>% 
        filter(qr.rank == max(qr.rank)) %>% 
        filter(is.effect | levels == max(levels)) %>% 
        pull(covariate) %>% 
        unique()
    #message('==============')
    # print(full.mod.matrix %>% t())
    # print(colinearities)
    # message(non.colinear.variables)
    # message(colinear.variables)
    #message('==============')
    if (effect.variable %in% colinear.variables) {
        non.colinear.variables
    } else {
        unique(c(non.colinear.variables, sample(colinear.variables, 1)))
    }
}

estimate_svs <- function(
    sample_metadata,
    counts_matrix,
    full_vars,
    reduced_vars=NULL,
    bind=TRUE,
    ...){
    # Run SVA to identify SVs and return SV matrix
    # reduced_vars=useful_covariates; full_vars=c(effect.variable, useful_covariates)
    # First convert variable lists to formula objs
    reduced_fml <- 
        {
            # set null model to include covariate(s) to remove
            if (is.null(reduced_vars) | (length(reduced_vars) == 0)) {
                "~ 1"
            } else {
                reduced_vars %>% 
                paste(collapse="+") %>% 
                sprintf('~%s', .)
            } 
        } %>% 
        formula() 
    # effects  + covaraites
    full_fml <- 
        full_vars %>% 
        paste(collapse="+") %>% 
        sprintf('~%s', .) %>% 
        formula() 
    # Run SVA 
    tryCatch(
        {
            svs <- 
                counts_matrix %>%
                as.matrix() %>%
                {.[rowSums(.) > 0, ]} %>% 
                svaseq(
                    # full model matrix
                    full_fml %>% model.matrix(sample_metadata),
                    # covariate only model matrix
                    reduced_fml %>% model.matrix(sample_metadata)
                ) %>%
                {.$sv} %>% 
                set_colnames(paste0('SV', 1:ncol(.)))
            if (bind) bind_cols(sample_metadata, svs) else svs
        },
        error=function(cond){
            message('Estimating SVs failed with following error:')
            message(conditionMessage(cond))
            if (bind) sample_metadata  else NA
        }
    )
}

run_deseq <- function(
    sample_metadata,
    counts_matrix,
    include.SVs,
    intercept,
    effect.variable,
    covariate.variables){
    # First check that covariates are useful (not uniform && not co-linear with effects)
    # message(paste(c(effect.variable, covariate.variables), collapse=','))
    useful_covariates <- 
        if (all(is.na(covariate.variables))) {
            list()
        } else { 
            covariate.variables
        }
    #     check_covariates_valid(
    #         sample_metadata,
    #         effect.variable,
    #         covariate.variables
    #     )
    message(glue('Covariates: {paste(useful_covariates, collapse=",")}'))
    # Estimate SVs if there are covariate to adjust for
    if (include.SVs) {
        sample_metadata <- 
            sample_metadata %>% 
            estimate_svs(
                counts_matrix,
                c(effect.variable, useful_covariates),
                useful_covariates,
            )
    }
    svs <- 
        sample_metadata %>%
        colnames() %>%
        {grep('^SV[0-9]', ., value=TRUE)}
    # Specify model formula (with intercept) i.e. ~ 1+ effect.variable + useful_covariates + SVs
    model.formula <- 
        c(
            intercept,
            effect.variable,
            useful_covariates,
            svs
        ) %>% 
        paste(collapse='+') %>% 
        sprintf('~ %s', .) %>%
        formula()
   message(""); print(model.formula)
    # Now run DESeq2 model and get results with all effects + useful covariates + SVs
    DESeqDataSetFromMatrix(
        countData=counts_matrix,
        colData=sample_metadata,
        design=model.formula
    ) %>%
    DESeq(
        test='Wald',
        parallel=TRUE
    )
}

get_contrast_results <- function(
    dds,
    gene_id_column,
    name=NULL,
    contrast=NULL,
    ...){
    # Format DESeq2 model results for specified model term (variable) into tidy table
    { 
        if (!is.null(contrast)) {
            DESeq2::results(dds, contrast=contrast)
        } else if (!is.null(name)) {
            DESeq2::results(dds, name=name)
        } else {
            stop('No contrast specified, add contrast or name argument')
        }
    } %>% 
    as.data.frame() %>%
    rownames_to_column(var=gene_id_column) %>%
    dplyr::as_tibble() %>%
    select(
        !!sym(gene_id_column),
        pvalue,
        log2FoldChange,
        lfcSE
    )
}

list_all_contrasts <- function(
    model.variables,
    sample_metadata,
    samples.list=NULL,
    sampleID.variable=NULL){
    if (!is.null(samples.list)) {
        sample_metadata <- 
            sample_metadata %>%
            filter((!!sym(sampleID.variable)) %in% samples.list)
    }
    # List all variables in the model
    tibble(model.variable=model.variables) %>%
    # Only keep categorical variables
    filter(!is.na(model.variable)) %>% 
    select(where(~ is.factor(.x) | is.character(.x))) %>% 
    # Get all pairwise combinations of levels for each variable i.e.
    # 1 row per contrast used to fetch DEG results
    group_by(model.variable) %>% 
    reframe(
        sample_metadata %>% 
        select(all_of(model.variable)) %>% 
        distinct() %>% 
        lapply(as.character) %>% 
        deframe() %>%
        combn(m=2) %>%
        t() %>%
        as_tibble() %>% 
        set_names(c('numerator', 'denominator')),
    ) %>% 
    mutate(contrast.name=glue('{model.variable}-{numerator}-{denominator}'))
}

get_all_contrast_results <- function(
    dds,
    model.variables,
    gene_id_column,
    ...){
    dds %>%
    resultsNames() %>%
    tibble(name=.) %>% 
    # Structure contrast info from names 
    rowwise() %>% 
    mutate(
        variable= 
            str_extract(
                name,
                paste(model.variables, collapse='|')
            ),
        comparison=
            ifelse(
                is.na(variable),
                NA,
                str_remove(name, glue('{variable}_'))
            ),
        variable=ifelse(is.na(variable), name, variable)
    ) %>%
    separate_wider_delim(
        comparison,
        delim='_vs_',
        names=
            c(
                'numerator', 
                'denominator'
            )
    ) %>% 
    mutate(
        deg.results=
            pmap(
                .l=.,
                get_contrast_results,
                gene_id_column=gene_id_column,
                dds=dds
            ) 
    )
}
###############
# Generate Bootstrapping
generate_bootstrap_sample_lists <- function(
    sample.metadata,
    effect.variable,
    n_bootstraps,
    chunk_size,
    n.samples.per.group,
    seed,
    sampleID.variable,
    ...){
    # sample.metadata=sample.metadata; effect.variable=EFFECT_VARAIBLE; sampleID.variable=SAMPLE_ID_VARIABLE; output_dir=BOOTSTRAP_LIST_DIR; n_bootstraps=N_BOOTSTRAPS; chunk_size=CHUNK_SIZE; n.samples.per.group=N_SAMPLES; seed=SEED
    set.seed(seed)
    # dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
    # Names of sample groups being compared
    groups <-
        sample.metadata[[effect.variable]] %>%
        unique() %>%
        as.character() %>%
        sort()
    if (length(groups) != 2) { stop('Comparison must have 2 groups only') }
    # comparison <- glue('{effect.variable}_{groups[[1]]}_vs_{groups[[2]]}')
    # Samples in group 1
    all_group1_samples <-
        sample.metadata %>%
        filter((!!sym(effect.variable)) == groups[1]) %>%
        pull(!!sym(sampleID.variable))
    # Samples in group 2
    all_group2_samples <-
        sample.metadata %>%
        filter((!!sym(effect.variable)) == groups[2]) %>%
        pull(!!sym(sampleID.variable))
    # Group bootstraps into chunks (files)
    n_chunks <- round(n_bootstraps / chunk_size)
    rep(
        1:n_chunks,
        chunk_size
    ) %>%
    paste0('chunk', .) %>% 
    tibble(chunk_id=.) %>%
    arrange(chunk_id) %>% 
    head(n=n_bootstraps) %>% 
    add_column(
        bs_id=paste0('bs', 1:n_bootstraps),
        comparison=glue('{effect.variable}_{groups[1]}_vs_{groups[2]}')
    ) %>%
    # for each bootstrap randomly pick samples for each group with replacement
    rowwise() %>%
    mutate(
        group1_samples=list(sample(all_group1_samples, size=n.samples.per.group)),
        group2_samples=list(sample(all_group2_samples, size=n.samples.per.group))
    ) %>% 
    mutate(across(ends_with('_samples'), \(x) paste(x, collapse=','))) %>% 
    ungroup()
}

generate_bootstraps <- function(
    model.specifications,
    sample.metadata){
    all.bootstrap.sample.lists <- 
        model.specifications %>% 
        mutate(
            bootstrap.sample.lists=
                purrr::pmap(
                    .,
                    generate_bootstrap_sample_lists,
                    sample.metadata=sample.metadata,
                    .progress=TRUE
                )
        ) %>% 
        unnest(bootstrap.sample.lists) %>% 
        # Now every row is 1 bootstrap, total n_bootstraps * n_models rows
        rowwise() %>%
        mutate(
            model.info=
                list(
                    glue('seed={seed}'),
                    glue('n.samples.per.group={n.samples.per.group}'),
                    glue('filter={filter.variable}'),
                    glue('CPM.cutoff={CPM.cutoff}'),
                    glue('min.frac.samples.expressed={min.frac.samples.expressed}'),
                    glue('intercept={intercept}'),
                    glue('include.SVs={include.SVs}'),
                    glue('effect={effect.variable}'),
                    glue('covariates={paste(covariate.variables, collapse=',')}')
                ) %>%
                list()
        ) %>% 
        ungroup() %>% 
        # collapse bootstraps into chunks
        nest(
            bootstrap.sample.lists=
                c(
                    bs_id,
                    comparison,
                    group1_samples,
                    group2_samples
                )
        ) %T>% 
        pmap(
            .l=.,
            .f=
                function(
                    list.dir, 
                    bootstrap.sample.lists,
                    chunk_id,
                    ...){
                    dir.create(
                        list.dir,
                        showWarnings=FALSE,
                        recursive=TRUE
                    )
                    write_tsv(
                        bootstrap.sample.lists,
                        file.path(list.dir, glue('{chunk_id}.tsv'))
                    )
                },
            .progress=TRUE
        )
    # Save conditions
    all.bootstrap.sample.lists %>% 
        distinct(
            info.filepath, 
            .keep_all=TRUE
        ) %>% 
        pmap(
            .l=.,
            .f=
                function(
                    model.info,
                    info.filepath,
                    model.filename, 
                    model.dir,
                    ...) {
                    sink(info.filepath)
                    paste(
                        'ARGUMENTS',
                        model.info,
                        sep='::'
                    ) %>% 
                    paste(collapse='\n') %>%
                    cat('\n')
                    print(session_info())
                    sink()
                    # copy model file for referencing
                    file.copy(
                        model.filename,
                        file.path(model.dir, 'this.model'),
                        overwrite=TRUE,
                        recursive=FALSE
                    )
                }
        )
}

deseq_pipeline <- function(
    sample_metadata,
    counts_matrix,
    sampleID.variable,
    filter.variable,
    CPM.cutoff,
    min.frac.samples.expressed,
    intercept,
    include.SVs,
    effect.variable,
    covariate.variables,
    gene_id_column,
    ...){ 
    # Filter low expression genes
    genes_to_keep <- 
        sample_metadata %>% 
        select(
            all_of(
                c(
                    sampleID.variable,
                    filter.variable
                )
           )
        ) %>%
        filter_genes_by_CPM(
            counts_matrix=counts_matrix,
            sampleID.variable=sampleID.variable,
            gene_id_column=gene_id_column,
            lib_sizes=NULL,
            CPM.cutoff=CPM.cutoff, 
            min.frac.samples.expressed=min.frac.samples.expressed 
        ) 
   message(glue('{length(genes_to_keep)} of {nrow(counts_matrix)} ({round(100 * length(genes_to_keep) / nrow(counts_matrix), 2)}%) genes retained after CPM filtering'))
    # Build DESeq2 model with specified parameters + samples on filtered genes
    run_deseq(
        sample_metadata=sample_metadata,
        counts_matrix=counts_matrix[genes_to_keep, ],
        include.SVs=include.SVs,
        intercept=intercept,
        effect.variable=effect.variable,
        covariate.variables=covariate.variables
    ) %>%
    # Fetch DEG results for all contrasts (named list of tables, 1 table per contrast)
    get_all_contrast_results(
        sample_metadata=sample_metadata,
        model.variables=c(effect.variable, covariate.variables),
        gene_id_column=gene_id_column
    )
}

run_pipeline_on_bootstrap <- function(
    sample.metadata,
    counts.matrix,
    sampleID.variable,
    samples.list,
    bs_id,
    results.dir,
    ...){
    # Only keep samples in this bootstrap
    message(     "|=========================================|")
    message(glue("| Starting pipeline for bootstrap {bs_id} |"))
    message(     "|=========================================|")
    # Subset metadata and counts to samples in this specific bootstrap
    sample_metadata <- 
        sample.metadata %>%
        filter((!!sym(sampleID.variable)) %in% samples.list) %>%
        mutate(across(where(is.factor), droplevels))
   message(glue('{nrow(sample_metadata)} samples'))
    counts_matrix <- 
        counts.matrix[, sample_metadata[[sampleID.variable]]]
   message(glue('{nrow(counts_matrix)} genes'))
    # Run entire pipeline to generate tables of DESeq2 results, 1 per contrast
    deseq_pipeline(
        sample_metadata=sample_metadata,
        counts_matrix=counts_matrix,
        sampleID.variable=sampleID.variable,
        ...
    ) %>%
    add_column(bs_id=bs_id, results.dir=results.dir)
}

run_pipeline_on_chunk <- function(
    bootstrap.chunk,
    chunk_id,
    model.parameters,
    sample.metadata,
    counts.matrix,
    sampleID.variable,
    gene_id_column,
    force_redo=FALSE,
    ...){
    # chunk_id=chunks.list$chunk_id[[1]]; bootstrap.chunk=chunks.list$bootstrap.chunk[[1]]; sampleID.variable=options$SAMPLE_ID_VAR; gene_id_column=options$GENE_ID_COLUMN; force_redo=options$FORCE_REDO; filter.variable=model.parameters$filter.variable[[1]]; CPM.cutoff=model.parameters$CPM.cutoff[[1]]; min.frac.samples.expressed=model.parameters$min.frac.samples.expressed[[1]]; intercept=model.parameters$intercept[[1]]; include.SVs=model.parameters$include.SVs[[1]]; effect.variable=model.parameters$effect.variable[[1]]; covariate.variables=model.parameters$covariate.variables[[1]]; 
    # samples.list=c(bootstrap.chunk$group1_samples[[3]],bootstrap.chunk$group2_samples[[3]]); 
    # check for cached results before computnig stuff
    results.files <- 
        list.files(
            model.parameters$results.dir[[1]],
            pattern=glue('{chunk_id}.tsv'),
            recursive=TRUE
        )
    if (length(results.files) > 0) {
        if (force_redo) {
            message('cached results exist, recomputing results anyways...')
            paste(results.files, sep='\n')
        } else {
            message('cached results exist, not recomputing')
            paste(results.files, sep='\n')
            return(NULL)
        }
    }
    # list all samples in each bootstrap
    message(tic())
    bootstrap.chunk %>% 
    rowwise() %>% 
    mutate(samples.list=list(c(group1_samples, group2_samples))) %>%
    # Add model parameters for all models specified
    cross_join(model.parameters) %>% 
    # build DESeq2 model + get all contrasts per bootstrap
    future_pmap(
        .l=.,
        # quietly(run_pipeline_on_bootstrap),
        run_pipeline_on_bootstrap,
        sample.metadata=sample.metadata,
        counts.matrix=counts.matrix,
        gene_id_column=gene_id_column,
        sampleID.variable=sampleID.variable,
        .progress=TRUE,
        .options=
            furrr_options(
                seed=model.parameters$seed[[1]],
                stdout=TRUE,
                chunk_size=5
             )
    ) %>%
    bind_rows() %>%
    # Combine bootstraps results into a single table, 1 column per bootstrap i.e.
    # produce 1 table  (PER contrast) combining all bootstraps in the chunk
    # | GeneName | statistic | bs1  | bs2  | bs3  | ... | bsN  |
    # |----------|-----------|------|------|------| ... |------|
    # | gene_A   | LFC       | 1.1  | 1.2  | 0.7  | ... | 1.1  |
    # | gene_A   | pvalue    | 0.1  | 0.7  | 0.01 | ... | 0.04 |
    # | gene_B   | LFC       | 2.1  | 2.2  | 1.7  | ... | 1.1  |
    # | gene_B   | pvalue    | 0.05 | 0.01 | 0.01 | ... | 0.4  |
    # ...
    rowwise() %>%
    mutate(contrast.name=glue('{variable}-{numerator}-{denominator}')) %>% 
    mutate(
        deg.results=
            deg.results %>%
            pivot_longer(
                cols=-all_of(c(gene_id_column)),
                names_to='stat',
                values_to=bs_id
            ) %>%
            list()
    ) %>% 
    group_by(contrast.name) %>% 
    summarize( 
        across(
            c(results.dir), # per contrast
            unique
        ),
        chunk.results=
            purrr::reduce(
                deg.results,
                full_join,
                by=c(gene_id_column, 'stat')
            ) %>%
            list()
    ) %>%
    # now save all combined bootstrap results in this chunk to 1 file per contrast
    mutate(
        results_file=
            file.path(
                results.dir,
                contrast.name,
                glue('{chunk_id}.tsv')
            )
    ) %>%
    pwalk(
        function(chunk.results, results_file, ...){
            dir.create(
                dirname(results_file),
                recursive=TRUE,
                showWarnings=FALSE
            )
            write_tsv(
                chunk.results,
                results_file,
            )
        }
    )
    message(toc())
}
###############
# Load Bootstraps
load_chunk <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(show_col_types=FALSE) %>% 
    pivot_longer(
        starts_with('bs'),
        names_to='bs_id',
        values_to='value'
    ) %>%
    pivot_wider(
        names_from=stat,
        values_from=value
    )
}

summarize_chunk <- function(
    filepath,
    adj.method='BH',
    ...){
    # filepath=file.path(BASE_DIR, 'results/bootstraps/SVA.only/DEG.results/Genotype-XDP-Control/chunk1.tsv'); adj.method='BH'
    filepath %>%
    load_chunk() %>% 
    # adjust pvalues within each bootstrap separately i.e. per ~20K genes
    group_by(bs_id) %>% 
    mutate(fdr=p.adjust(pvalue, method='BH')) %>% 
    ungroup() %>%
    # mutate(direction=ifelse(log2FoldChange > 0, 'up', 'down')) %>% 
    # Now compute summary stats+frequencies for each gene across all bootstraps in this chunk
    group_by(EnsemblID) %>% 
    summarize(
        # LFC stats
        meanLFC=mean(log2FoldChange, na.rm=TRUE),
        medianLFC=median(log2FoldChange, na.rm=TRUE),
        medianLFCSE=median(lfcSE, na.rm=TRUE),
        varLFC=var(log2FoldChange, na.rm=TRUE),
        # quantile stats
        max.fdr=max(fdr),
        min.fdr=min(fdr),
        # count number of bootstraps with various qualities
        n.attempts.either=n(),
        n.any.either=sum(!is.na(log2FoldChange)),
        n.any.up=sum(log2FoldChange > 0),
        n.any.down=sum(log2FoldChange < 0),
        # nominal counts
        n.nominal.either=sum(pvalue < 0.05),
        n.nominal.up=    sum(pvalue < 0.05 & log2FoldChange > 0),
        n.nominal.down=  sum(pvalue < 0.05 & log2FoldChange < 0),
        # fdr counts
        n.fdr.either=sum(fdr < 0.1),
        n.fdr.up=    sum(fdr < 0.1 & log2FoldChange > 0),
        n.fdr.down=  sum(fdr < 0.1 & log2FoldChange < 0)
    )

}

summarize_all_chunks <- function(
    chunk.files,
    ...){
    # chunks.dir=file.path(BASE_DIR, 'results/bootstraps/SVA.only/DEG.results')
    chunk.files %>% 
    tibble(filepath=.) %>%
    mutate(chunk_id=str_remove(basename(filepath), '.tsv')) %>%
    mutate(
        chunk.summarizes=
            future_pmap(
                 .l=filepath,
                 .f=summarize_chunk,
                 ...,
                 .progress=FALSE
            )
    ) %>%
    unnest(chunk.summarizes) %>%
    # Now agg. statistics across chunks
    group_by(EnsemblID) %>%
    summarize(
        meanLFC=mean(meanLFC, na.rm=TRUE),
        medianLFC=mean(medianLFC, na.rm=TRUE),
        medianLFCSE=median(medianLFCSE, na.rm=TRUE),
        # pooled variance across equally sized subgroups
        varLFC=mean(varLFC, na.rm=TRUE) + var(meanLFC, na.rm=TRUE),
        # var.lfc=var(var.lfc),
        # min over all mins
        across(
            starts_with('min.'),
            min
        ),
        # max over all maxs
        across(
            starts_with('max.'),
            max
        ),
        # sum counts
        across(
            starts_with('n.'),
            sum
        )
    ) %>% 
    # now calculate frequencies from each count variable
    mutate(
        across(
            starts_with('n.'),
            ~ .x / n.any.either,
            .names="freq.{.col}"
        )
    ) %>% 
    rename_with(~ str_replace(.x, '^freq.n.', 'freq.')) %>%
    mutate(freq.attempts.either=1)
# bootstrap.summary.df %>% select(model.name, contrast, n.attempts.either, n.any.either)
}

summarize_all_contrasts <- function(
    model.parameters,
    model.names=NULL,
    contrasts=NULL,
    force_redo=FALSE){
    # specify input, output files for each model + contrast
    model.parameters %>%
    # 1 row per model per contrast (i.e. per DEG results)
    rowwise() %>% 
    mutate(contrast=list(list.dirs(results.dir, recursive=FALSE, full.names=FALSE))) %>% 
    unnest(contrast) %>%
    # Filter unspecified models+contrasts
    {
        if (!is.null(model.names)) {
            filter(., model.name %in% model.names)
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(contrasts)) {
            filter(., contrast %in% contrasts)
        } else {
            .
        }
    } %>% 
    rowwise() %>% 
    mutate(
        chunk.files=
            file.path(results.dir, contrast) %>% 
            list.files(
                full.names=TRUE,
                pattern='chunk[0-9]+.tsv',
                recursive=FALSE
            ) %>%
            list(),
        results_file=
            file.path(
                results.dir, 
                glue('{contrast}.bootstrap.summary.tsv')
            )
    ) %>% 
    select(
        -c(
            model.filename,
            ends_with('.filepath'),
            ends_with('.dir')
        )
    ) %>% 
    # separate_wider_delim(
    #     contrast,
    #     delim='-',
    #     names=
    #         c(
    #             'variable',
    #             'numerator',
    #             'denominator'
    #         ),
    #     cols_remove=FALSE
    # ) %>% 
        # {.} %>% select(contrast, chunk.files, results_file)
    # now load + compute summary statistics on each chunk of bootstraps
    ungroup() %>% 
    mutate(
        results=
            pmap(
                .l=.,
                check_cached_results,
                force_redo=force_redo,
                return_data=TRUE,
                results_fnc=summarize_all_chunks,
                .progress=TRUE
            )
    ) %>% 
    nest(
        model.params=
            -c(
                model.name,
                contrast,
                results
            )
    ) %>% 
    unnest(results)
}
###############
# Analysis
subsample_bootstraps <- function(
    results.dir,
    model_name,
    contrast_name,
    n_subsamples=5,
    seed=9){
    set.seed(seed)
    results.dir %>% 
    list.files(
        full.names=TRUE,
        pattern='chunk[0-9]+.tsv',
        recursive=FALSE
    ) %>%
    sample(size=n_subsamples) %>% 
    tibble(filepath=.) %>%
    mutate(chunk_id=str_remove(basename(filepath), '.tsv')) %>%
    add_column(
        model.name=model_name,
        contrast=contrast_name
    ) %>% 
    rowwise() %>% 
    mutate(chunk.results=list(summarize_all_chunks(filepath))) %>%
    ungroup() %>% 
    unnest(chunk.results) %>% 
    select(
        c(
            model.name,
            contrast,
            EnsemblID,
            chunk_id,
            starts_with('freq.')
        )
    ) %>%
    pivot_longer(
        -c(
            model.name,
            contrast,
            chunk_id,
            EnsemblID
        ),
        names_to='stat',
        values_to='frequency'
    ) %>%
    separate_wider_delim(
        'stat',
        delim='.',
        names=
            c(
                'tmp',
                'significance',
                'direction'
            )
    ) %>%  
    select(-c(tmp))
}

join_full_and_bootstrap_results <- function(
    bootstrap.summary.df,
    full.model.results){
    bootstrap.summary.df %>%
    select(
        c(
            model.name,
            contrast,
            EnsemblID,
            starts_with('freq.')
        )
    ) %>%
    pivot_longer(
        -c(
            model.name,
            contrast,
            EnsemblID
        ),
        names_to='stat',
        values_to='frequency'
    ) %>%
    separate_wider_delim(
        'stat',
        delim='.',
        names=
            c(
                'tmp',
                'significance',
                'direction'
            )
    ) %>%  
    select(-c(tmp)) %>% 
    filter(
        significance != 'attempts',
        !(significance == 'any' & direction == 'either')
    ) %>% 
    inner_join(
        full.model.results %>%
        select(
            model.name,
            contrast,
            EnsemblID,
            log2FoldChange,
            pvalue,
            fdr
        ),
        by=
            join_by(
                model.name,
                contrast,
                EnsemblID
            )
    )
}

