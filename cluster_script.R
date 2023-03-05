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
cfg <- didehpc::didehpc_config(cluster = "wpia-hn",
                               template = 'AllNodes',
                               cores = 29)
# cfg <- didehpc::didehpc_config(cluster = "big",
#                                template = '32Core',
#                                cores = 29)
obj <- didehpc::queue_didehpc(ctx, config = cfg)


## Long run of the fits task
test_L.6.1 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.2 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.3 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.4 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))

test_L.6.5 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.6 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.7 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.8 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.9 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                               parameters = list(region = "london",
                                                                 short_run = TRUE,
                                                                 deterministic = TRUE,
                                                                 N = 1000,
                                                                 method = 'f1',
                                                                 start = 0,
                                                                 step = 0.01,
                                                                 ne = 25,
                                                                 L = 6),
                                               use_draft = "newer"))


test_L.6.10 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                                parameters = list(region = "london",
                                                                  short_run = TRUE,
                                                                  deterministic = TRUE,
                                                                  N = 1000,
                                                                  method = 'f1',
                                                                  start = 0,
                                                                  step = 0.01,
                                                                  ne = 25,
                                                                  L = 6),
                                                use_draft = "newer"))
test_L.6.1$status()
test_L.6.2$status()
test_L.6.3$status()
test_L.6.4$status()
test_L.6.5$status()
test_L.6.6$status()
test_L.6.7$status()
test_L.6.8$status()
test_L.6.9$status()
test_L.6.10$status()




test <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient',
                                         parameters = list(region = "london",
                                                           short_run = TRUE,
                                                           deterministic = TRUE,
                                                           N = 1000,
                                                           method = 'f1',
                                                           start = 0.36,
                                                           step = 0.04,
                                                           ne = 6,
                                                           L = 1),
                                         use_draft = "newer"))

test_trans.10 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                         parameters = list(region = "london",
                                                           short_run = TRUE,
                                                           deterministic = TRUE,
                                                           N = 1000,
                                                           method = 'f1',
                                                           start = 0.37,
                                                           step = 0.03,
                                                           ne = 10,
                                                           L = 1,
                                                           scale = 10),
                                         use_draft = "newer"))

test_trans.50 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                  parameters = list(region = "london",
                                                                    short_run = TRUE,
                                                                    deterministic = TRUE,
                                                                    N = 1000,
                                                                    method = 'f1',
                                                                    start = 0.37,
                                                                    step = 0.03,
                                                                    ne = 10,
                                                                    L = 1,
                                                                    scale = 50),
                                                  use_draft = "newer"))

test_trans.100 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                  parameters = list(region = "london",
                                                                    short_run = TRUE,
                                                                    deterministic = TRUE,
                                                                    N = 1000,
                                                                    method = 'f1',
                                                                    start = 0.37,
                                                                    step = 0.03,
                                                                    ne = 10,
                                                                    L = 1,
                                                                    scale = 100),
                                                  use_draft = "newer"))

test_trans.200 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                  parameters = list(region = "london",
                                                                    short_run = TRUE,
                                                                    deterministic = TRUE,
                                                                    N = 1000,
                                                                    method = 'f1',
                                                                    start = 0.37,
                                                                    step = 0.03,
                                                                    ne = 10,
                                                                    L = 1,
                                                                    scale = 200),
                                                  use_draft = "newer"))

test_trans.300 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                  parameters = list(region = "london",
                                                                    short_run = TRUE,
                                                                    deterministic = TRUE,
                                                                    N = 1000,
                                                                    method = 'f1',
                                                                    start = 0.37,
                                                                    step = 0.03,
                                                                    ne = 10,
                                                                    L = 1,
                                                                    scale = 300),
                                                  use_draft = "newer"))

test_trans.500 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                   parameters = list(region = "london",
                                                                     short_run = TRUE,
                                                                     deterministic = TRUE,
                                                                     N = 1000,
                                                                     method = 'f1',
                                                                     start = 0.37,
                                                                     step = 0.03,
                                                                     ne = 10,
                                                                     L = 1,
                                                                     scale = 500),
                                                   use_draft = "newer"))

test_trans.500$status()
test_trans.300$status()
test_trans.200$status()
test_trans.100$status()
test_trans.50$status()
test_trans.10$status()

test_trans.1000.1 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                    parameters = list(region = "london",
                                                                      short_run = TRUE,
                                                                      deterministic = TRUE,
                                                                      N = 1000,
                                                                      method = 'f1',
                                                                      start = 0.27,
                                                                      step = 0.03,
                                                                      ne = 16,
                                                                      L = 1,
                                                                      scale = 1000),
                                                    use_draft = "newer"))

test_trans.1000.2 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                      parameters = list(region = "london",
                                                                        short_run = TRUE,
                                                                        deterministic = TRUE,
                                                                        N = 1000,
                                                                        method = 'f1',
                                                                        start = 0.27,
                                                                        step = 0.03,
                                                                        ne = 16,
                                                                        L = 1,
                                                                        scale = 1000),
                                                      use_draft = "newer"))

test_trans.1000.3 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                      parameters = list(region = "london",
                                                                        short_run = TRUE,
                                                                        deterministic = TRUE,
                                                                        N = 1000,
                                                                        method = 'f1',
                                                                        start = 0.27,
                                                                        step = 0.03,
                                                                        ne = 16,
                                                                        L = 1,
                                                                        scale = 1000),
                                                      use_draft = "newer"))

test_trans.1000.4 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                      parameters = list(region = "london",
                                                                        short_run = TRUE,
                                                                        deterministic = TRUE,
                                                                        N = 1000,
                                                                        method = 'f1',
                                                                        start = 0.27,
                                                                        step = 0.03,
                                                                        ne = 16,
                                                                        L = 1,
                                                                        scale = 1000),
                                                      use_draft = "newer"))

test_trans.1000.5 <- obj$enqueue(orderly::orderly_run('vaccine_delay_fits_gradient_transform',
                                                      parameters = list(region = "london",
                                                                        short_run = TRUE,
                                                                        deterministic = TRUE,
                                                                        N = 1000,
                                                                        method = 'f1',
                                                                        start = 0.27,
                                                                        step = 0.03,
                                                                        ne = 16,
                                                                        L = 1,
                                                                        scale = 1000),
                                                      use_draft = "newer"))


test_trans.1000.5$status()
test_trans.1000.4$status()
test_trans.1000.3$status()
test_trans.1000.2$status()
test_trans.1000.1$status()