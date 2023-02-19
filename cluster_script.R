options(didehpc.cluster = "fi--didemrchnb",
        didehpc.username = "wl5818")

setwd(orderly::orderly_config()$root)

packages <- c("sircovid", "orderly", "mcstate", "dust", "spimalot", "dplyr","MASS","readr","ggplot2")

src <- conan::conan_sources(NULL,
                            repos = c("https://ncov-ic.github.io/drat",
                                      "https://raphaels1.r-universe.dev"))
ctx <- context::context_save("contexts",
                             packages = packages,
                             package_sources = src)
cfg <- didehpc::didehpc_config(cluster = "big",
                               template = '32Core',
                               cores = 32)
obj <- didehpc::queue_didehpc(ctx, config = cfg)


## Long run of the fits task
fit <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                        parameters = list(region = "london",
                                                          short_run = TRUE,
                                                          deterministic = TRUE,
                                                          N = 1000,
                                                          method = 'f1',
                                                          start = 0.04,
                                                          step = 0.01,
                                                          ne = 15,
                                                          L = 2),
                                        use_draft = "newer")
)
fit$result()

