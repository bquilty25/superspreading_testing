# Genuine Prose Changes Versus manuscript.md

This note compares the legacy manuscript at `docs/manuscript.md` with the current manuscript rendered from `docs/manuscript.qmd` to Markdown. It excludes passages where the wording is effectively unchanged and only the evaluated numbers differ.

## 1. Contact-duration interval label

Old (`docs/manuscript.md`):

> Contact durations in CoMix differed significantly between household and non-household contacts, with a median duration of 30 minutes (IQR: 5, 180 minutes) for out-of-household contacts compared to 480 minutes (8 hours) (95% CI: 180, 1080 minutes) for household contacts (Figure S1).

New (current rendered manuscript):

> Contact durations in CoMix differed significantly between household and non-household contacts, with a median duration of 30 minutes (IQR: 5, 180 minutes) for out-of-household contacts compared to 480 minutes (8 hours) (IQR: 180, 1080 minutes) for household contacts (Figure S1).

Assessment:

This is a genuine editorial correction: the household interval label changes from `95% CI` to `IQR`.

## 2. Discussion contact-summary sentence

Old (`docs/manuscript.md`):

> Daily contact rates became both lower on average (from ~12 per day to ~6 per day) and more overdispersed during the pandemic in the UK in 2020 as some individuals maintained high contact rates (e.g., essential workers) while others had their number of daily contacts reduced to zero (those working from home).

New (current rendered manuscript):

> Daily contact rates became both lower on average, from 11.5 per day to 4.7 per day, and more overdispersed during the pandemic in the UK in 2020 as some individuals maintained high contact rates, such as essential workers, while others had their number of daily contacts reduced to zero, such as those working from home.

Assessment:

This is mostly a move from approximate to exact values, but the sentence was also lightly rewritten to remove the parenthetical shorthand.

## 3. Discussion contact-heterogeneity paragraph

Old (`docs/manuscript.md`):

> Contact heterogeneity was found to contribute more to superspreading than heterogeneity in viral load, as the number of contacts an individual makes places an upper limit on the number of people they can infect, and, despite some individuals having much higher infectious potential than others, most infected individuals pass through a high viral load period (~68% with a peak culture probability >0.6, ~52% with a peak culture probability >0.8).. Viral load heterogeneity does still contribute significantly to heterogeneity in transmission since viral load varies by orders of magnitude over time within individuals and between individuals (Figure 5). This indicates that superspreading is more a result of between-person contact variability and within-person infectiousness variability than between-person infectiousness variability. Hence if reductions in contact rates can be targeted specifically at individuals ahead or during their window of high infectivity (e.g., by encouraging contacts of cases to undergo regular self-testing and self-isolation upon a positive rapid test result 22) then this may result in reductions in transmission while minimising the burden of quarantine.

New (current rendered manuscript):

> Contact heterogeneity was found to contribute more to superspreading than heterogeneity in viral load, as the number of contacts an individual makes places an upper limit on the number of people they can infect, and, despite some individuals having much higher infectious potential than others, most infected individuals pass through a high viral load period, with 68% having a peak culture probability greater than 0.6 and 53% greater than 0.8. Viral load heterogeneity does still contribute significantly to heterogeneity in transmission since viral load varies by orders of magnitude over time within individuals and between individuals (Figure 5). This indicates that superspreading is more a result of between-person contact variability and within-person infectiousness variability than between-person infectiousness variability. Hence if reductions in contact rates can be targeted specifically at individuals ahead of or during their window of high infectivity, for example by encouraging contacts of cases to undergo regular self-testing and self-isolation upon a positive rapid test result [@quilty2021], then this may result in reductions in transmission while minimising the burden of quarantine.

Assessment:

This is a real prose clean-up. The current manuscript fixes punctuation, repairs the awkward parenthetical phrasing, and smooths the transition into the intervention sentence.

## 4. Discussion testing paragraph

Old (`docs/manuscript.md`):

> Like others 6,26,30, we hypothesised that lateral flow testing, by detecting individuals with high viral loads when they were most infectious, would reduce transmission through reducing the potential for superspreading. This manifested as a decrease in the proportion infecting over 10 others and a substantial increase in the proportion infecting zero others. This, perhaps counter-intuitively, resulted in a decrease in k in some cases, as the relative increase in those infecting zero others exceeded the decrease in those infecting many others (here, over 10 others). Hence, assessment of superspreading solely via the metric of the overdispersion parameter k may conceal changes in both the upper and lower tail of the secondary infection distribution. Both regular testing and pre-event testing were effective in reducing R given high enough frequency or a low enough event size threshold, respectively, as long as uptake or adherence was high. Testing had the highest relative impact on transmission when contact rates were high (e.g., at pre-pandemic levels) as there were more potentially preventable exposures, meaning rapid testing could reduce R below the growth threshold of 1 while otherwise maintaining relatively normal contact rates. In contrast, testing during lockdown would have less impact as R was already below 1. Having everyone in the population test every 3 days would bring R below 1 and be approximately equivalent in terms of impact on secondary infections to having people test only before events of minimum size 10 for pre-pandemic contact levels. This indicates that rapid testing could be an effective, minimally disruptive intervention to reduce transmission if uptake/adherence could be maximised through incentivising use.

New (current rendered manuscript):

> Like others [@ke2022; @hart2023; @middleton2024], we hypothesised that lateral flow testing, by detecting individuals with high viral loads when they were most infectious, would reduce transmission through reducing the potential for superspreading. This manifested as a decrease in the proportion infecting over 10 others and a substantial increase in the proportion infecting zero others. This, perhaps counter-intuitively, resulted in a decrease in $k$ in some cases, as the relative increase in those infecting zero others exceeded the decrease in those infecting many others, here defined as over 10 others. Hence, assessment of superspreading solely via the metric of the overdispersion parameter $k$ may conceal changes in both the upper and lower tail of the secondary infection distribution. Both regular testing and pre-event testing were effective in reducing $R$ given high enough frequency or a low enough event size threshold, respectively, as long as uptake or adherence was high. Testing had the highest relative impact on transmission when contact rates were high, for example at pre-pandemic levels, as there were more potentially preventable exposures, meaning rapid testing could reduce $R$ below the growth threshold of 1 while otherwise maintaining relatively normal contact rates. In contrast, testing during lockdown would have less impact as $R$ was already below 1. Having everyone in the population test every 3 days would bring $R$ below 1 and be approximately equivalent in terms of impact on secondary infections to having people test only before events of minimum size 10 for pre-pandemic contact levels. This indicates that rapid testing could be an effective, minimally disruptive intervention to reduce transmission if uptake or adherence could be maximised through incentivising use.

Assessment:

This paragraph is mostly the same argument, but it has been copy-edited: citation style is normalised, `k` and `R` are formatted consistently, and several phrases are made more explicit.

## 5. Final concluding sentence

Old (`docs/manuscript.md`):

> Our results suggest superspreading for SARS-CoV-2 can be best explained as a random sample from the tail of the contact and shedding distribution: it occurs when an infected individual makes a high number of contacts during a highly infectious period lasting approximately 2 days on average, with the majority of infected individuals being sufficiently infectious (~68% with a peak culture probability over 0.6, ~52% over 0.8) least one day to be capable of causing a superspreading event -given they make a high number of contacts.

New (current rendered manuscript):

> Our results suggest superspreading for SARS-CoV-2 can be best explained as a random sample from the tail of the contact and shedding distribution: it occurs when an infected individual makes a high number of contacts during a highly infectious period lasting approximately 2 days on average, with the majority of infected individuals being sufficiently infectious, 68% with a peak culture probability over 0.6 and 53% over 0.8, for at least one day to be capable of causing a superspreading event, given they make a high number of contacts.

Assessment:

This is the clearest genuine prose repair. The old sentence was grammatically broken; the current manuscript restores it to a readable and logically complete conclusion.