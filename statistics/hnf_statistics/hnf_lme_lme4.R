# Called by hnf_lme_octave on GNU Octave releases without fitlmematrix.
# Each model includes a fixed and random intercept. The independent random
# slope terms fit a diagonal covariance for the numeric, pre-expanded design.
fit_maps <- function(folder) {
    if (!requireNamespace("lme4", quietly = TRUE))
        stop("Install the R package lme4 or use Octave statistics fitlmematrix.")
    input <- file(file.path(folder, "data.bin"), "rb")
    on.exit(close(input), add = TRUE)
    read_int <- function(n) readBin(input, integer(), n, size = 4, endian = "little")
    read_double <- function(n) readBin(input, numeric(), n, size = 8, endian = "little")
    dims <- read_int(6)
    n <- dims[1]
    voxels <- dims[2]
    p <- dims[3]
    slopes <- read_int(dims[4])
    reml <- as.logical(dims[5])
    diagonal <- as.logical(dims[6])
    subject <- factor(read_double(n))
    design <- as.data.frame(matrix(read_double(n*p), n, p))
    terms <- paste0("x", seq_len(p))
    names(design) <- terms
    design$subject <- subject
    # Eagerly reading all responses is unnecessary. Columns are consecutive
    # in the input, so only one response vector is resident in R at a time.
    separator <- if (diagonal && length(slopes)) " || " else " | "
    random <- paste(c("1", terms[slopes]), collapse = " + ")
    model <- reformulate(c(terms, paste0("(", random, separator, "subject)")), "y")
    beta <- tstat <- matrix(NA_real_, p, voxels)
    singular <- 0L
    for (v in seq_len(voxels)) {
        design$y <- read_double(n)
        if (length(design$y) != n) stop("Truncated response data")
        fit <- suppressMessages(lme4::lmer(model, design, REML = reml))
        singular <- singular + lme4::isSingular(fit)
        coef <- coef(summary(fit))
        beta[,v] <- coef[terms, "Estimate"]
        tstat[,v] <- coef[terms, "t value"]
    }
    if (singular) message(singular, " of ", voxels, " mixed fits on the covariance boundary")
    output <- file(file.path(folder, "fit.bin"), "wb")
    on.exit(close(output), add = TRUE)
    writeBin(as.double(beta), output, size = 8, endian = "little")
    writeBin(as.double(tstat), output, size = 8, endian = "little")
}
fit_maps(Sys.getenv("NEUROFREQ_LME_DIR"))
