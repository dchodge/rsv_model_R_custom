prop_mat <- read.table(here("inst", "extdata", "demo", "prop_par.txt"))[[1]]
prop_risk <- read.table(here("inst", "extdata", "demo", "proprisk.txt"))
pVHR <- prop_risk[1:25, ]
pHR <- prop_risk[26:50, ]
pLR <- prop_risk[51:75, ]

cnt_matrix_p <- read.table(here("inst", "extdata", "cnt_mat", "cntPA.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
cnt_matrix_p_h <- read.table(here("inst", "extdata", "cnt_mat", "cntPAH.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
pwp_p <- read.table(here("inst", "extdata", "cnt_mat", "cntPpwp.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
pwn_p <- read.table(here("inst", "extdata", "cnt_mat", "cntPpwn.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
nwp_p <- read.table(here("inst", "extdata", "cnt_mat", "cntPnwp.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
nwn_p <- read.table(here("inst", "extdata", "cnt_mat", "cntPnwn.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)


cnt_matrix_c <- read.table(here("inst", "extdata", "cnt_mat", "cntCA.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
cnt_matrix_c_h <- read.table(here("inst", "extdata", "cnt_mat", "cntCAH.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
pwp_c <- read.table(here("inst", "extdata", "cnt_mat", "cntCpwp.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
pwn_c <- read.table(here("inst", "extdata", "cnt_mat", "cntCpwn.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
nwp_c <- read.table(here("inst", "extdata", "cnt_mat", "cntCnwp.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)
nwn_c <- read.table(here("inst", "extdata", "cnt_mat", "cntCnwn.txt"))[[1]] %>% matrix(25, 25, byrow = TRUE)

data_inter_uk <- list(prop_mat = prop_mat, 
    pVHR = pVHR,
    pHR = pHR,
    pLR = pLR,
    cnt_matrix_p = cnt_matrix_p,
    cnt_matrix_p_h = cnt_matrix_p_h,
    pwp_p = pwp_p,
    pwn_p = pwn_p,
    nwp_p = nwp_p,
    nwn_p = nwn_p,
    cnt_matrix_c = cnt_matrix_c,
    cnt_matrix_c_h = cnt_matrix_c_h,
    pwp_c = pwp_c,
    pwn_c = pwn_c,
    nwp_c = nwp_c,
    nwn_c = nwn_c
    )

save(data_inter_uk, file = here("data", "inter_data_uk.RData"))