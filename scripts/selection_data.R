library(dplyr)

con <- src_sqlite("raw_data/Experiment_Data.db", create=FALSE)

obswin <- tbl(con, "ObsWindow") %>% collect(n = Inf)
tdata <- tbl(con, "TrialData") %>% collect(n = Inf)

obsdat <- tdata %>% inner_join(obswin, "RespID") %>%
    mutate(CritSel=ObjSel=="unme",
           PGd=Cond=="PG",
           CGd=Cond=="CG",
           PG=PGd-mean(PGd), CG=CGd-mean(CGd)) %>%
    select(SubjID, ItemID, RespID, PG, CG, Cond, RT=rtArtOn,
           CritSel, ObjSel)

cond_lookup <- data_frame(Cond=c("CG", "PG", "NG"),
                          Condition=c("Common Ground",
                              "Privileged Ground",
                              "Control"))

obsel_lookup <- data_frame(ObjSel=c("unme", "comp", "negt"),
                           Alternative=c("Unnamed",
                               "Critical",
                               "Negated"))

inner_join(obsdat, cond_lookup) %>%
    inner_join(obsel_lookup) %>%
    group_by(Condition, Alternative) %>%
    summarize(Y=n()) %>%
    inner_join(obsdat %>% inner_join(cond_lookup) %>%
                   inner_join(obsel_lookup) %>%
                   group_by(Condition) %>%
                   summarize(N=n())) %>%
    mutate(p=round(Y/N, 3))

library(lme4)

obsdat.glmer <-
    glmer(CritSel ~ PG + CG + (PG + CG | SubjID) + (PG + CG | ItemID),
          data=obsdat, family=binomial(link=logit),
          control=glmerControl(optimizer="bobyqa"))

obsdat.glmer.noME <-
    glmer(CritSel ~ (PG + CG | SubjID) + (PG + CG | ItemID),
         data=obsdat, family=binomial(link=logit),
         control=glmerControl(optimizer="bobyqa"))

obsdat.glmer.noCG <-
    glmer(CritSel ~ PG + (PG + CG | SubjID) + (PG + CG | ItemID),
          data=obsdat, family=binomial(link=logit),
          control=glmerControl(optimizer="bobyqa"))

obsdat.glmer.noPG <-
    glmer(CritSel ~ CG + (PG + CG | SubjID) + (PG + CG | ItemID),
          data=obsdat, family=binomial(link=logit),
          control=glmerControl(optimizer="bobyqa"))

summary(obsdat.glmer)
anova(obsdat.glmer, obsdat.glmer.noME)
anova(obsdat.glmer.noCG, obsdat.glmer)
anova(obsdat.glmer.noPG, obsdat.glmer)
