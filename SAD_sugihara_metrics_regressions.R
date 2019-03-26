
#############
require(tidyverse)
require(MuMIn)
require(lmerTest)

#############

metrics.results <- readr::read_delim(file = "./results/sugihara_model_metrics_complete.csv",delim = ";")

resource.dist.levels <- unique(metrics.results$resource.distribution.level)
richness.levels <- unique(metrics.results$richness.level)
connectance.levels <- unique(metrics.results$connectance.level)
apportionment.levels <- unique(metrics.results$apportionment.level)
replicates <- max(metrics.results$replicate)
trophic.guilds <- unique(metrics.results$trophic.guild)


metrics.results$resource.distribution.level <- factor(metrics.results$resource.distribution.level,levels = c("uniform","skewed"))
metrics.results$richness.level <- factor(metrics.results$richness.level)
metrics.results$connectance.level <- factor(metrics.results$connectance.level)
metrics.results$apportionment.level <- recode(metrics.results$apportionment.level,"random.assortment" = "random assortment",
                                              "random.fraction" = "random fraction",
                                              "dominance.decay" = "dominance decay",
                                              "dominance.preemption" = "dominance preemption")
# reorder levels
metrics.results$apportionment.level <- factor(metrics.results$apportionment.level,levels = c("random assortment","dominance decay","random fraction","dominance preemption"))
# metrics.results$trophic.guild <- factor(metrics.results$trophic.guild, levels = c("primary Producers","herbivores","omnivores","carnivores"))
metrics.results$trophic.guild <- factor(metrics.results$trophic.guild,levels = c("resources","basal","intermediate","top"))

#####

metrics.model.data <- droplevels(subset(metrics.results,trophic.guild != "resources" & resource.distribution.level == "skewed"))

# null.model <- lmer(hill.evenness ~ vulnerability.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
# tl.model <- lmer(hill.evenness ~ vulnerability.level + trophic.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
# 
# anova(null.model,tl.model)

evenness.sugihara <- lmerTest::lmer(hill.evenness ~ apportionment.level + connectance.level + trophic.guild + (1 | replicate), data = metrics.model.data, REML = FALSE)
evenness.fixed <- lm(hill.evenness ~ apportionment.level + connectance.level + trophic.guild,data = metrics.model.data)
evenness.fixed.int <- lm(hill.evenness ~ apportionment.level+connectance.level+trophic.guild+apportionment.level*connectance.level,data = metrics.model.data)

# skewness.sugihara <- lmerTest::lmer(skewness ~  richness.level + connectance.level + trophic.guild + (1 | replicate), data = metrics.model.data, REML = FALSE)

MuMIn::r.squaredGLMM(evenness.sugihara)
# MuMIn::r.squaredGLMM(skewness.sugihara)

summary(evenness.sugihara)

# in this case, the linear model is virtually identical to the mixed effects one, so I can stick with the simpler
# for calculating the proportion of variance explained by each factor, this is simply the percentage of the sum of squares

model.ssq <- anova(evenness.fixed)$"Sum Sq"
model.ssq/sum(model.ssq)*100

# alternatively, check the r2 of nested models including the random factor
# both approaches are almost identical, so it's ok

evenness.m2 <- lmerTest::lmer(hill.evenness ~ apportionment.level + connectance.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
MuMIn::r.squaredGLMM(evenness.m2)

evenness.m3 <- lmerTest::lmer(hill.evenness ~ apportionment.level + (1 | replicate), data = metrics.model.data, REML = FALSE)
MuMIn::r.squaredGLMM(evenness.m3)

