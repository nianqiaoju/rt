## rm(list = ls())
set.seed(1)
## setwd() if necessary

## load libraries
library(ggplot2)
library(dplyr)
library(devtools)
devtools::install_github("qingyuanzhao/bets.covid19")
library(bets.covid19)
library(EpiEstim)

## color blind palette for ggplot
cbPalette <- c("#0072B2","#D55E00", "#E69F00","#999999", "#56B4E9", "#009E73", "#F0E442",   "#CC79A7")
## red, blue, orange, grey,  light blue, green, yellow, pink

## read in the raw counts and process the data
df <- read.csv("ExtractCounts.csv")
df <- df[-which(is.na(df$Count)),]
df$Raw <- NULL
df$Date <- as.Date(df$Date, "%d-%b")
df$Date[df$Date > as.Date("2020-04-30")] <- df$Date[df$Date > as.Date("2020-04-30")] - 366
print(head(df))

## Visualize the symptom curve
ggplot(df) + aes(x = Date, y = Count) + geom_bar(stat = "identity") + theme_bw()

## first reproduce the Rt plot using Cori et al.
## see https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
## it uses all the data but only displays result in Janurary
It <- df$Count
## it uses 5-day moving average
## starts on December 21th, which is the 14th entry
t_start <- seq(14, length(It) - 5)
t_end <- t_start + 5
## "The daily number of reported COVID-19 cases and
## the serial interval (mean, 7.5 days [SD, 3.4 days]; constant across periods),
## derived from a previous epidemiological survey of the first 425 cases in Wuhan,[8]
## were used to estimate Rt and its 95% credible interval on each day via a 5-day moving average."
cori_res_parametric_si <- estimate_R(It,
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 7.5,
                                  std_si = 3.4,
                                  t_start = t_start,
                                  t_end = t_end)))
# print(head(res_parametric_si$R))
# plot(res_parametric_si, legend = FALSE)
## only plot results since Janurary
num_days <- length(cori_res_parametric_si$R$t_end)
results <- data.frame(method = 'Cori et al.', date = integer(num_days), meanR = double(num_days), lbd = double(num_days), ubd = double(num_days))
## t_start starts at 14, which correspond to the date 21-Dec-2019
results$date <- df$Date[cori_res_parametric_si$R$t_end]
results$meanR <- cori_res_parametric_si$R$`Mean(R)`
results$lbd <- cori_res_parametric_si$R$`Quantile.0.025(R)`
results$ubd <- cori_res_parametric_si$R$`Quantile.0.975(R)`
head(results %>% subset(date >= as.Date('2020-01-01')))

wt_res_parametric_si <- wallinga_teunis(It,
                                        method = 'parametric_si',
                                        config = list(
                                            mean_si = 7.5,
                                            std_si = 3.4,
                                            t_start = t_start,
                                            t_end = t_end,
                                            n_sim = 100))
print(wt_res_parametric_si$R)
wt_results <- data.frame(method = 'Wallinga-Teunis (smoothing)',
           date = df$Date[wt_res_parametric_si$R$t_end],
           meanR = wt_res_parametric_si$R$`Mean(R)`,
           lbd = wt_res_parametric_si$R$`Quantile.0.025(R)`,
           ubd = wt_res_parametric_si$R$`Quantile.0.975(R)`)
wt_results_january <- subset(wt_results, date >= 32)
results <- rbind(results, wt_results)

## No smoothing (t_start = t_end)
wt_res_parametric_si <- wallinga_teunis(It,
                                        method = 'parametric_si',
                                        config = list(
                                            mean_si = 7.5,
                                            std_si = 3.4,
                                            t_start = t_end,
                                            t_end = t_end,
                                            n_sim = 100))
print(wt_res_parametric_si$R)
wt_results <- data.frame(method = 'Wallinga-Teunis (no smoothing)',
           date = df$Date[wt_res_parametric_si$R$t_end],
           meanR = wt_res_parametric_si$R$`Mean(R)`,
           lbd = wt_res_parametric_si$R$`Quantile.0.025(R)`,
           ubd = wt_res_parametric_si$R$`Quantile.0.975(R)`)
wt_results_january <- subset(wt_results, date >= 32)
results <- rbind(results, wt_results)


## break_days <- sapply(c('1-Jan','10-Jan','23-Jan', '2-Feb','17-Feb'), date.process) + as.Date("2019-11-30")
## print(break_days)
## p1 <- ggplot(data = subset(results, date >= as.Date('2020-01-01')), aes(x = date, y = meanR, group = method))
## p1 <- p1 + geom_ribbon(aes(ymin = lbd, ymax = ubd, fill = method), alpha = 0.3) + geom_line( aes(colour = method))
## p1 <- p1 + geom_vline(xintercept = break_days)
## p1 <- p1 + geom_hline(yintercept = 1)
## p1 <- p1 + scale_y_continuous(minor_breaks = seq(0,5,0.5)) + scale_x_date(date_breaks = "5 days", date_minor_breaks = "1 day", date_labels = "%b %d")
## p1 <- p1 + scale_fill_manual(values = cbPalette) + scale_colour_manual(values=cbPalette)
## p1 <- p1+ ylab('Rt') + xlab('Date') + theme_bw(base_size = 18) + theme(legend.position = "bottom")
## p1

## ggsave('reproduce_figure4.pdf', width = 12, height = 8)


############################################################################################
############################## Estimate Rt using simualted It ##############################
############################################################################################

rt_start <- as.Date("2020-01-01")
rt_end <- as.Date("2020-02-29")

simulate.It <- function(St, mean_ip = 5.6, std_ip = 4.1) {

    alpha <- (mean_ip / std_ip)^2
    beta <- mean_ip / std_ip^2

    It <- data.frame(Date = seq(St$Date[1] - round(qgamma(0.9999, alpha, beta)), max(St$Date), 1),
                     Count = 0)

    for (i in 1:nrow(St)) {
        if (St$Count[i] > 0) {
            ip <- floor(rgamma(St$Count[i], alpha, beta) + 0.5)
            ip <- table(ip)
            for (j in 1:length(ip)) {
                day <- St$Date[i] - as.numeric(names(ip)[j])
                It$Count[It$Date == day] <- It$Count[It$Date == day] + ip[j]
            }
        }
    }

    It
}

St <- df

simulate.Rt <- function(sim) {
    It <- simulate.It(St)

    t_start = seq(which(It$Date == as.Date("2019-12-21")), nrow(It) - 5)
    t_end = t_start + 5

    res <- estimate_R(It$Count,
                      method="parametric_si",
                      config = make_config(list(
                          mean_si = 7.5,
                          std_si = 3.4,
                          t_start = t_start,
                          t_end = t_end)))

    output <- data.frame(sim = sim,
                         Date = seq(as.Date("2019-12-21") + 5, It$Date[nrow(It)], 1),
                         Rt = sample_posterior_R(res, length(t_start), 1:length(t_start)))

    output
}

output <- parallel::mclapply(1:1000, simulate.Rt, mc.cores = 3)
output <- do.call(rbind, output)

## break_days <- sapply(c('1-Jan','10-Jan','23-Jan', '2-Feb','17-Feb'), date.process) + as.Date("2019-11-30")
## p2 <- ggplot(subset(output, Date >= as.Date("2020-01-01")))+
##     aes(x = Date, y = Rt) + geom_line(aes(group = sim), alpha = 0.3) +
##     geom_line(stat = "summary", fun.y = "mean", col = "red") +
##     theme_bw()
## p2 <- p2 + geom_vline(xintercept = break_days)
## p2 <- p2 + geom_hline(yintercept = 1)
## p2 <- p2 + scale_y_continuous(minor_breaks = seq(0,5,0.5))
## p2 <- p2 + scale_x_date(date_breaks = "5 days", date_minor_breaks = "1 day", date_labels = "%b %d")
## p2 <- p2 + ylab('Rt') + xlab('Date')
## p2

## ggsave('rt_use_simulated_it.pdf', width = 12, height = 8)

## ## compare the WT curve with the curved genereted by (Cori + inferred It)
## wt_curve <- results[which(results$method == 'Wallinga-Teunis (smoothing)'),]$meanR
## wt_curve
## cori_i_curve <- data.frame(output %>% group_by(Date) %>% summarise(meanR = mean(Rt)))
## cori_i_curve
## plot(cori_i_curve, type = 'l', col = 'red', ylim = c(0,4))
## lines(cori_i_curve$Date, wt_curve)

## organize the output data into results format and compare with cori and wt
infer_i_result <- output %>% group_by(Date) %>% summarise(method = 'Cori et al. (back-calculated incidence)',meanR = mean(Rt), lbd = quantile(Rt, 0.025), ubd = quantile(Rt, 0.975))
head(infer_i_result)
names(infer_i_result) <- c('date','method','meanR','lbd','ubd')
results <- rbind(results, infer_i_result)
print(tail(results))

results$method <- factor(results$method, levels = c("Cori et al.", "Wallinga-Teunis (smoothing)", "Cori et al. (back-calculated incidence)", "Wallinga-Teunis (no smoothing)"))

p3 <- ggplot(data = results, aes(x = date, y = meanR, group = method))
#p3 <- ggplot(data = subset(results, date >= as.Date('2020-01-01')), aes(x = date, y = meanR, group = method))
p3 <- p3 + geom_ribbon(aes(ymin = lbd, ymax = ubd, fill = method), alpha = 0.3) + geom_line( aes(colour = method, linetype = method))
p3 <- p3 + geom_vline(xintercept = break_days)
p3 <- p3 + geom_hline(yintercept = 1)
p3 <- p3 + scale_y_continuous(minor_breaks = seq(0,5,0.5)) + scale_x_date(date_breaks = "5 days", date_minor_breaks = "1 day", date_labels = "%b %d")
p3 <- p3 + scale_fill_manual(values = cbPalette[c(1,2,1,2)]) + scale_colour_manual(values=cbPalette[c(1,2,1,2)]) + scale_linetype_manual(values = c("dashed", "dashed", "solid", "solid"))
p3 <- p3 + ylab('R(t)') + xlab('Date') + theme_bw(base_size = 18) + theme(legend.position = "bottom") + labs(color = "Method", fill = "Method", linetype = "Method") + guides(fill=guide_legend(nrow=2))
p3

ggsave('rt_four_methods.pdf', width = 12, height = 8)
ggsave('rt_four_methods.png', width = 12, height = 8)
