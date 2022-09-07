# libraries we need
libs <- c(
    "reshape2", "tidyverse",
    "quanteda", "rstan", "parallel"
)

# install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
    install.packages(libs[!installed_libs])
}

# load libraries
invisible(lapply(libs, library, character.only = T))

#   Dim 2: Intervention (a.k.a. AvP)
#   Juraj Medzihorsky
#   2017-02-28
#   misc
rstan_options(auto_write = TRUE)
options(mc.cores = 3)
options(stringsAsFactors = F)


#   load the data with the speeches
ukr_unsc <- read.csv("ukr_unsc.csv") %>% select(-Agenda)
df <- melt(ukr_unsc, id.vars = c("Source", "Date"))
names(df) <- c("source", "date", "country", "speech")
d <- df[, 1:4] %>%
    group_by(source) %>%
    mutate(meeting = cur_group_id()) %>%
    ungroup() %>%
    select(meeting, date, country, speech)

head(d)

#   process the data
#   load and process the dictionary
dictionary_text <- readLines("dictionary.txt", warn = F)

dict <- sapply(
    dictionary_text, function(x) strsplit(x, ": ")[[1]][-1],
    USE.NAMES = F
)
dict <- strsplit(dict, ", ")[-5]
names(dict) <- c("proc", "act", "conf", "crim")

#   create the DFM
mat_1 <- dfm(d$speech, dictionary = dictionary(dict), stem = TRUE)

#   work copy
pm <- as.data.frame(as.matrix(mat_1))

#   a new dataframe, for now all countries
a <- d[, 1:4]
a$wp1 <- pm$proc
a$wp2 <- pm$conf
a$wn1 <- pm$act
a$wn2 <- pm$crim
a$n1 <- pm$proc + pm$act
a$n2 <- pm$conf + pm$crim


#   long form, only dim CvS here
l <- a[a$n2 > 0, ]

#   find countries that spoke at > 1 meeting
countries_to_take <- names(table(l$country)[table(l$country) > 1])

l <- l[l$country %in% countries_to_take, ]

#   recode the meetings into integers by overall chronological order
meetings <- unique(l$meeting)
meetings <- meetings[order(as.numeric(meetings))]
l$meet_fact <- match(l$meeting, meetings)

#   first split the data by country
s <- split(l, f = as.factor(l$country))

#   a function to recode within countries
aux_relative <- function(x) {
    x <- x[order(x$meet_fact), ]
    x$relative <- 1:nrow(x)
    x$relat_max <- nrow(x)
    x$dist <-
        c(
            999^2,
            sapply(
                2:nrow(x), function(i) x$meet_fact[i] - x$meet_fact[i - 1]
            )
        )
    x$dist[x$dist == 0] <- 1 # account for double meetings
    x$sqrt_dist <- sqrt(x$dist)
    return(x)
}

#   apply the function and merge the data
ll <- l
ll <- do.call(rbind, lapply(s, aux_relative))

#   countries into factors
countries <- unique(ll$country)
countries <- countries[order(countries)]

ll$coun_fact <- match(ll$country, countries)

#   prepare the data for stan
ldx <- list(
    y = ll$wp2,
    nw = ll$n2,
    meet_fact = ll$meet_fact,
    meet_relative = ll$relative,
    country = ll$coun_fact,
    n_ms = nrow(ll),
    n_c = length(countries),
    n_m = length(meetings)
)
country_frame <- matrix(nrow = ldx$n_c, ncol = 2)
for (tc in 1:ldx$n_c) {
    country_frame[tc, ] <- range(which(ldx$country == tc))
}
ldx$country_frame <- country_frame
ldx$omega <- l$sqrt_dist


model_code_1 <- "
    data {
        int<lower=1> n_ms;
        int<lower=1> n_c;
        int<lower=1> n_m;
        int<lower=0> y[n_ms];
        int<lower=1> nw[n_ms];
        int<lower=1> meet_fact[n_ms];
        int<lower=1> meet_relative[n_ms];
        int<lower=1> country[n_ms];
        int country_frame[n_c, 2];
        real<lower=1> omega[n_ms];
    }
    parameters {
        vector[n_m] gamma_star;
        real<lower=0> sigma_gamma;
        vector<lower=0>[n_c] sigma_theta;
        real<lower=0> tau;
        vector<lower=-pi()/2, upper=pi()/2>[n_ms] theta_step_unif;
        vector[n_ms] theta_norm;
    }
    transformed parameters {
        vector[n_ms] theta;
        for (i in 1:n_ms) {
            if (meet_relative[i]==1) {
                theta[i] = theta_norm[i];
            } else {
                theta[i] = theta[i-1] + tan(theta_step_unif[i])*omega[i]*sigma_theta[country[i]]*tau;
            }
        }
    }
    model {
        y ~ binomial_logit(nw, gamma_star[meet_fact]*sigma_gamma + theta);
        theta_step_unif ~ uniform(-pi()/2, pi()/2);
        theta_norm ~ normal(0, 1);
        sigma_theta ~ normal(0, 1);
        tau ~ normal(0, 1);
        gamma_star ~ normal(0, 1);
        sigma_gamma ~ normal(0, 1);
    }
    generated quantities {
        real country_mean[n_c];
        real country_sd[n_c];
        vector[n_m] gamma;
        vector[n_c] daily;
        for (tc in 1:n_c) {
            country_mean[tc] = mean(theta[country_frame[tc, 1]:country_frame[tc, 2]]);
            country_sd[tc] = sd(theta[country_frame[tc, 1]:country_frame[tc, 2]]);
        }
        for (tm in 1:n_m)
            gamma[tm] = gamma_star[tm] * sigma_gamma;
        for (tc in 1:n_c)
            daily[tc] = sigma_theta[tc]*tau;
    }"


r_test <- stan(
    model_code = model_code_1, data = ldx, chains = 3, iter = 1e1 * 4,
    control = list(
        adapt_delta = 0.99,
        max_treedepth = 15
    )
)

print(r_test, pars = c("daily[1]", "tau", "gamma[1]", "sigma_gamma", "lp__"))

pairs(r_test, pars = c("daily[1]", "tau", "gamma[1]", "sigma_gamma", "lp__"), las = 1)

st <- summary(r_test)
summary(st$summary[, "Rhat"])

#   rename objects and save them
dim_cvs <- r_test
data_cvs <- ldx
save(list = c("countries", "data_cvs", "dim_cvs"), file = "sim_cvs.RData")
