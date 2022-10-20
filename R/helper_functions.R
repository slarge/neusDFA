
sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}


dfa_mod <- function(dat, m, R, cov_v = NA, just_testing = TRUE, data_wide = TRUE, safe = FALSE, ...)  {


  if(just_testing) {
    # cntl.list <- list(minit = 200, maxit = 5200, allow.degen = TRUE, conv.test.slope.tol = 0.5)
    cntl.list = list(minit = 200,
                     maxit = 40000, allow.degen = TRUE,
                     abstol = 0.0001, conv.test.slope.tol = 0.5,
                     safe = safe)
  }
  if(!just_testing){
    cntl.list <- list(minit = 200, maxit = 60000, allow.degen = FALSE,
                      abstol = 0.0001, conv.test.slope.tol = 0.1,
                      safe = safe)
  }
  if(!data_wide){
    dat <- dat %>%
      tidyr::pivot_wider(names_from = year_season, values_from = est) %>%
      tibble::column_to_rownames(var = "spp") %>%
      as.matrix
  }

  mod_list <- list(m = m, R = R)

  if(length(cov_v) == 1 & all(is.na(cov_v))){
    m1 <- MARSS(dat,
                model = mod_list,
                form = "dfa",
                z.score = FALSE,
                control = cntl.list)
  }
  if(length(cov_v) > 1 & all(!is.na(cov_v))){

    # cov_v <- as.vector(cov_v)


    m1 <- MARSS(dat,
                model = mod_list,
                form = "dfa",
                z.score = FALSE,
                covariates = cov_v,
                control = cntl.list)
  }
  return(m1)
}



# dfa_mod <- function(dat, m, R, covariate, just_testing = TRUE, ...)  {
#
#
#   if(just_testing) {
#     cntl.list = list(minit = 200, maxit = 3000, allow.degen = TRUE)
#   }
#   if(!just_testing){
#     cntl.list <- list(minit = 200, maxit = 60000, allow.degen = FALSE,
#                       abstol = 0.0001, conv.test.slope.tol = 0.1)
#   }
#   ## control parameters if just testing convergence tolerance
#
#   #
#   # cntl.list <- list(minit = 200, maxit = 5200, allow.degen = TRUE, conv.test.slope.tol = 0.01)
#   mod_list <- list(m = m, R = R)
#
#   if(covariate == "none"){
#     m1 <- MARSS(temp_wide,
#                 model = mod_list,
#                 form = "dfa",
#                 z.score = FALSE,
#                 control = cntl.list)
#   }
#   if(covariate == "season_f"){
#     # mod_list <- list(m = m, R = R)
#     m1 <- MARSS(temp_wide,
#                 model = mod_list,
#                 form = "dfa",
#                 z.score = FALSE,
#                 covariates = season_f,
#                 control = cntl.list)
#   }
#   return(m1)
# }






Z_maker <- function(n, m, nsites = 1) {
  if(!n%%nsites == 0)
    stop("nsites is not a multiple of n")
  n <- n/nsites

  # Set up default Z
  Z <- matrix(list(), nrow = n, ncol = m)
  # insert row (i) & col (j) indices
  for (i in seq(n)) {
    Z[i, ] <- paste("Z", i, seq(m), sep = "")
  }
  # set correct i,j values in Z to numeric 0
  if (m > 1) {
    for (i in 1:(m - 1)) {
      Z[i, (i + 1):m] <- 0
    }
  }

  return(do.call(rbind, replicate(nsites, Z, simplify = FALSE)))
}



#' Reload a VAST model
#'
#' \code{reload_model} allows a user to save a fitted model, reload it in a new
#'      R terminal, and then relink the DLLs so that it functions as expected.
#'
#' @inheritParams make_model
#' @param x Output from \code{\link{fit_model}}, potentially with DLLs not linked
#' @param check_gradient Whether to check the gradients of the reloaded model
#'
#' @return Output from \code{\link{fit_model}} with DLLs relinked
#'
#' @examples
#' \dontrun{
#' # Run model
#' fit = fit_model( ... )
#' saveRDS( object=fit, file="path_and_name.rds" )
#'
#' # Reload and relink
#' fit_new = readRDS( file="path_and_name.rds" )
#' fit_new = reload_model( x = fit_new )
#' }
#'
#' @export
reload_model <-
  function( x,
            check_gradient = TRUE,
            CompileDir = system.file("executables",package = "VAST"),
            Version = x$settings$Version,
            framework = x$input_args$model_args_input$framework,
            Obj = x$tmb_list$Obj ){

    # Load old one
    if( is.null(framework) ){
      Version_framework = Version
    }else{
      Version_framework = paste0( Version, "_", framework )
    }
    origwd = getwd()
    on.exit( setwd(origwd), add=TRUE )
    setwd(CompileDir)
    dyn.load( TMB::dynlib(Version_framework) ) # random=Random,

    # Retape
    Obj$retape()

    # Ensure that last.par and last.par.best are right
    Obj$fn(x$parameter_estimates$par)

    # Check gradient
    if( check_gradient==TRUE ){
      Gr = Obj$gr(x$parameter_estimates$par)
      if(max(abs(Gr))>1){
        warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": does not seem converged")
      }else if(max(abs(Gr))>0.01){
        warning("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": might not be converged")
      }else{
        message("Maximum absolute gradient of ", signif(max(abs(Gr)),3), ": No evidence of non-convergence")
      }
    }

    return(x)
  }


Z_maker <- function(n, m, nsites = 1) {
  if(!n%%nsites == 0)
    stop("nsites is not a multiple of n")
  n <- n/nsites

  # Set up default Z
  Z <- matrix(list(), nrow = n, ncol = m)
  # insert row (i) & col (j) indices
  for (i in seq(n)) {
    Z[i, ] <- paste("Z", i, seq(m), sep = "")
  }
  # set correct i,j values in Z to numeric 0
  if (m > 1) {
    for (i in 1:(m - 1)) {
      Z[i, (i + 1):m] <- 0
    }
  }

  return(do.call(rbind, replicate(nsites, Z, simplify = FALSE)))
}



dfa_plot <- function(object, EPU = NULL) {

  marss_obj <- MARSSparamCIs(object)

  name_code <- rownames(object$marss$data)
  m <- nrow(object$states.se)

  if(!is.null(EPU)){
    # title_name <- switch(EPU,
    #                      GOM = "Gulf of Maine",
    #                      MAB = "Mid-Atlantic Bight",
    #                      SS = "Scotian Shelf",
    #                      GB = "Georges Bank")
    title_name = EPU
  }
  if(is.null(EPU)){
    title_name = ""
  }

  # the rotation matrix for the Z
  z <- coef(marss_obj, type = "Z")
  H.inv <- varimax(z)$rotmat

  # Get the Z, upZ, lowZ
  z.low <- coef(marss_obj, type = "Z", what="par.lowCI")
  z.up <- coef(marss_obj, type = "Z", what="par.upCI")

  z.rot <- z %*% H.inv
  z.rot.up <- z.up %*% H.inv
  z.rot.low <- z.low %*% H.inv


  factor_df <- data.frame(name_code = rep(name_code, m),
                          trend = rep(paste0("Trend ", 1:m), each = length(name_code)),
                          Z = as.vector(z),
                          Zup = as.vector(z.up),
                          Zlow = as.vector(z.low)) %>%
    mutate(spp = gsub("_",  " ", name_code),
           spp = ifelse(grepl("america?|atlant?|acadia?", spp),
                        stringr::str_to_sentence(spp),
                        spp),
           shape_id = 19, #ifelse(Zup < 0 | Zlow > 0, 19, 1),
           # fill_id = ifelse(Zup < 0 | Zlow > 0, "black", "grey70"),
           fill_id = ifelse(trend == "Trend 1", "#264CFF", "#FF420E"),
           color_id = fill_id,
           alpha_id = ifelse(Zup < 0 | Zlow > 0, .9, .8),
    ) %>%
    arrange(desc(spp)) %>%
    mutate(spp = factor(spp, levels = unique(spp)))



  fplot <-  ggplot(data = factor_df,
                   aes(x = spp,
                       ymin = Zlow,
                       ymax = Zup,
                       y = Z)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = 0.2, color = "grey60", linetype = "dashed") +
    geom_hline(yintercept = -0.2, color = "grey60", linetype = "dashed") +
    geom_pointrange(aes(shape = shape_id, fill = fill_id, color = color_id, alpha = alpha_id), show.legend = FALSE)+ #, position = position_dodge(width = 0.5)) +
    # facet_wrap(~trend) +
    coord_flip() +
    scale_shape_identity() +
    scale_fill_identity() +
    scale_color_identity() +
    labs(title = title_name,
         x = "", y = "factor loadings") +
    theme_minimal()

  trend_df <- tsSmooth(object, type = "xtT", interval = "confidence",
                       level = .95) %>%
    mutate(
      trend = gsub(pattern = "^X", "Trend ", .rownames),
      year = rep(as.numeric(colnames(object$marss$data)), m),
      fill_id = ifelse(trend == "Trend 1", "#264CFF", "#FF420E"),
      color_id = fill_id)


  tplot <- ggplot(data = trend_df,
                  aes(x = year,
                      ymin = .conf.low,
                      ymax = .conf.up,
                      y = .estimate)) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_line(aes(color = color_id)) +
    geom_ribbon(aes(fill = fill_id), alpha = 0.6, show.legend = FALSE)+ #, position = position_dodge(width = 0.5)) +
    # facet_wrap(~ trend, ncol = m) +
    labs(x = "", y = "Estimate") +
    scale_fill_identity() +
    scale_color_identity() +
    theme_minimal()

  tfplot <- fplot/tplot + plot_layout(heights = c(3, 1))
  #

  return(tfplot)

}