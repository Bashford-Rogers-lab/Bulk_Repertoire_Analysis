## Function to compare how the Correlation Matrix compares to Cor.Matrix generated using different methods 
## Allows you to chose the best method for your data
## Lauren Overend
## Jan 2023
## Lauren.overend@oriel.ox.ac.uk

#day1 <- readRDS('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Correlation_between_measures_SUBSAMPLED_DAY1ONLY_BCR_PRODUCTIVE.rds')
#day3 <- readRDS('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Correlation_between_measures_SUBSAMPLED_DAY3ONLY_BCR_PRODUCTIVE.rds')
#day5 <-readRDS('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Correlation_between_measures_SUBSAMPLED_DAY5ONLY_BCR_PRODUCTIVE.rds')
#allday <-readRDS('/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR/Summary/Correlation_between_measures_SUBSAMPLED_BCR_PRODUCTIVE.rds')
#outputdir <- '/gpfs2/well/immune-rep/shared/MISEQ/SEPSIS_COMPLETE/BCR'

compare_methods <- function(day1, day3, day5, allday, outputdir){
			
			allx <- c()
			#------------------------------------------------------
			#### CHANGE IN CORRELATION using DAY1 R 
			## Lets look at change in correlation - should be very small unless the values are deviating 
			Changeday35 <-day3- day1
			Changeday15 <-day5- day1
			
			mat_feature_corr_long_day3 <- data.frame(Changeday35)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "DCorr"
			mat_feature_corr_long_day3$DAY <- "3"
			mat_feature_corr_long_day3$METHOD <- "DAY1R"
			
			mat_feature_corr_long_day5 <- data.frame(Changeday15)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "DCorr"
			mat_feature_corr_long_day5$DAY <- "5"
			mat_feature_corr_long_day5$METHOD <- "DAY1R"
			
			## BIND
			allx <- rbind(allx, mat_feature_corr_long_day3, mat_feature_corr_long_day5)
			
			#------------------------------------------------------
			#### CHANGE IN CORRELATION using DAY3 R 
			## Lets look at change in correlation - should be very small unless the values are deviating 
			Changeday35 <-day1- day3
			Changeday15 <-day5- day3
			
			mat_feature_corr_long_day3 <- data.frame(Changeday35)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "DCorr"
			mat_feature_corr_long_day3$DAY <- "1"
			mat_feature_corr_long_day3$METHOD <- "DAY3R"
			
			mat_feature_corr_long_day5 <- data.frame(Changeday15)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "DCorr"
			mat_feature_corr_long_day5$DAY <- "5"
			mat_feature_corr_long_day5$METHOD <- "DAY3R"
			
			## BIND
			allx <- rbind(allx, mat_feature_corr_long_day3, mat_feature_corr_long_day5)
			
			#------------------------------------------------------
			#### CHANGE IN CORRELATION using DAY5 R 
			## Lets look at change in correlation - should be very small unless the values are deviating 
			Changeday35 <-day1- day5
			Changeday15 <-day3- day5
			
			mat_feature_corr_long_day3 <- data.frame(Changeday35)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "DCorr"
			mat_feature_corr_long_day3$DAY <- "1"
			mat_feature_corr_long_day3$METHOD <- "DAY5R"
			
			mat_feature_corr_long_day5 <- data.frame(Changeday15)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "DCorr"
			mat_feature_corr_long_day5$DAY <- "3"
			mat_feature_corr_long_day5$METHOD <- "DAY5R"
			
			## BIND
			allx <- rbind(allx, mat_feature_corr_long_day3, mat_feature_corr_long_day5)
			
			
			#------------------------------------------------------
			#####################################################################################
			### CHANGE IN CORRELATION using AVERAGE R 
			avgcorr <- (day1 +  day3+day5/3)
			Changeday13 <-day1- avgcorr
			Changeday35 <-day3- avgcorr
			Changeday15 <-day1- avgcorr
			
			mat_feature_corr_long_day1 <- data.frame(Changeday13)
			mat_feature_corr_long_day1$Chain <- rownames(mat_feature_corr_long_day1)
			mat_feature_corr_long_day1 <- mat_feature_corr_long_day1 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day1)[3] <- "DCorr"
			mat_feature_corr_long_day1$DAY <- "1"
			mat_feature_corr_long_day1$METHOD <- "AVGR"
			
			mat_feature_corr_long_day3 <- data.frame(Changeday35)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "DCorr"
			mat_feature_corr_long_day3$DAY <- "3"
			mat_feature_corr_long_day3$METHOD <- "AVGR"
			
			mat_feature_corr_long_day5 <- data.frame(Changeday15)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "DCorr"
			mat_feature_corr_long_day5$DAY <- "5"
			mat_feature_corr_long_day5$METHOD <- "AVGR"
			
			## RBIND
			allx <- rbind(allx, mat_feature_corr_long_day1, mat_feature_corr_long_day3, mat_feature_corr_long_day5)
			
			#####################################################################################
			### CHANGE IN CORRELATION using ALLDAYS R 
			Changeday13 <-day1- allday
			Changeday35 <-day3- allday
			Changeday15 <-day1- allday
			
			mat_feature_corr_long_day1 <- data.frame(Changeday13)
			mat_feature_corr_long_day1$Chain <- rownames(mat_feature_corr_long_day1)
			mat_feature_corr_long_day1 <- mat_feature_corr_long_day1 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day1)[3] <- "DCorr"
			mat_feature_corr_long_day1$DAY <- "1"
			mat_feature_corr_long_day1$METHOD <- "ALLR"

			mat_feature_corr_long_day3 <- data.frame(Changeday35)
			mat_feature_corr_long_day3$Chain <- rownames(mat_feature_corr_long_day3)
			mat_feature_corr_long_day3 <- mat_feature_corr_long_day3 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day3)[3] <- "DCorr"
			mat_feature_corr_long_day3$DAY <- "3"
			mat_feature_corr_long_day3$METHOD <- "ALLR"
			
			mat_feature_corr_long_day5 <- data.frame(Changeday15)
			mat_feature_corr_long_day5$Chain <- rownames(mat_feature_corr_long_day5)
			mat_feature_corr_long_day5 <- mat_feature_corr_long_day5 %>% gather(chain2, value, -c(Chain))
			colnames(mat_feature_corr_long_day5)[3] <- "DCorr"
			mat_feature_corr_long_day5$DAY <- "5"
			mat_feature_corr_long_day5$METHOD <- "ALLR"
			
			allx <- rbind(allx, mat_feature_corr_long_day1, mat_feature_corr_long_day3, mat_feature_corr_long_day5)

			#####################################################################################

			### lets summarise 
			## Generate some summary statistics describing the distribution
			stats <- allx  %>%group_by(METHOD, DAY) %>%dplyr::summarize(Mean = mean(DCorr, na.rm=TRUE), SD = sd(DCorr, na.rm=TRUE), Median= median(DCorr, na.rm=TRUE))
			stats <- data.frame(stats)
			stats1 <-  gather(stats, Statistic, value, Mean:Median)
			write.table(stats1, paste0(outputdir, "/Summary/R_Methods_Compare_Stats.txt"))
			## Plot 
			pdf(paste0(outputdir, "/Plots/R_Methods_Compare.pdf"), width=10, height=10)
			p1 <- ggplot(allx, aes(x=DCorr, fill=DAY, group=DAY), alpha=0.5)+geom_histogram(bins=50) +theme_classic()+ylab("Count")+xlab("Difference in R")+guides(fill="none")+geom_vline(xintercept=0.5, col="blue")+geom_vline(xintercept=-0.5, col="blue")+facet_grid(rows=vars(DAY), cols=vars(METHOD))
			p2 <- ggplot(allx, aes(y=DCorr, x=DAY, fill=DAY), alpha=0.5)+geom_violin() +theme_classic()+ylab("Difference in R")+xlab("Day used for generating Pearson Correlation Matrix")+guides(fill="none")+geom_hline(yintercept=0.5, col="blue")+geom_hline(yintercept=-0.5, col="blue")+facet_grid(cols=vars(METHOD))
			p3 <- ggplot(stats1, aes(y=value, x=METHOD, fill=DAY, colour=DAY))+geom_point(alpha=0.7) +theme_classic()+ylab("Value")+xlab("Difference in R")+facet_grid(cols=vars(Statistic))+
			scale_x_discrete(guide = guide_axis(n.dodge=2))+geom_hline(yintercept=0, col="blue")
			plot(plot_grid(p1, p2, p3, labels="AUTO", align="h", axis="lbt", ncol=1))
			dev.off()
			}
			
			
			