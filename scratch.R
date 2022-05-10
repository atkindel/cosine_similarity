# Get log-rank frequency ratio and length of difference vector
# Use 1% of Jeopardy clue token pairs
jtok_dists %>%
  slice_sample(prop=0.01) %>%
  rowwise() %>%
  mutate(rn_f = sum(str_detect(j_cluetext$clue_text, token), str_detect(j_cluetext$answer, token)),
         rn2_f = sum(str_detect(j_cluetext$clue_text, token2), str_detect(j_cluetext$answer, token2)),
         d_norm = norm(glove_d[rn,] - glove_d[rn2,], type="2")) %>%
  mutate(rnr = log(rn / rn2),
         rnr_b = log(rn * rn2),
         rnfr = log(rn_f / rn2_f)) ->
  test
 
test %>%
  ggplot(aes(x=cos_sim, y=d_norm, color=rnfr)) +
  geom_point()





# They share a basis
# test <- gb_bl_basis$u %*% t(bw_bl_basis$u) %*% bw_bl_basis$u %*% t(gb_bl_basis$u)
# test_basis <- svd(test)
# sum(test_basis$d)
diff_ev <- rbind(data.frame(x=1:length(gb_bl_basis$d), y=gb_bl_basis$d, group="Good words"),
                 data.frame(x=1:length(bw_bl_basis$d), y=bw_bl_basis$d, group="Bad words"))
diff_ev %>%
  ggplot(aes(x, y, color=group)) +
  geom_point() +
  geom_line()

d_prop_diff <- (gb_bl_basis$d / sum(gb_bl_basis$d)) / (bw_bl_basis$d / sum(bw_bl_basis$d))
data.frame(x=1:length(gb_bl_basis$d),
           y=(gb_bl_basis$d / sum(gb_bl_basis$d)) / (bw_bl_basis$d / sum(bw_bl_basis$d))) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  labs(x="Eigenvector (exemplar)")

min(d_prop_diff) / max(d_prop_diff)


rbind(data.frame(x=1:length(wwg_basis$d),
                 y=(wwg_basis$d / sum(wwg_basis$d)) / (wbg_basis$d / sum(wbg_basis$d)),
                 group="Names"),
      data.frame(x=1:length(good_basis$d),
                 y=(good_basis$d / sum(good_basis$d)) / (bad_basis$d / sum(bad_basis$d)),
                 group="Attributes"),
      data.frame(x=1:length(gb_bl_basis$d),
                 y=(gb_bl_basis$d / sum(gb_bl_basis$d)) / (bw_bl_basis$d / sum(bw_bl_basis$d)),
                 group="Euclidean diff.")) ->
  eigenratios

eigenratios %>%
  ggplot(aes(x, y, color=group)) +
  geom_hline(yintercept=1, linetype="dashed") + 
  geom_point() +
  labs(x="Principal component",
       y="Eigenratio",
       color="Subspace")

eigenratios %>%
  pivot_wider(names_from = group, values_from=y) ->
  ertest



The predicted relationship can be shown in word embeddings. As an example, I randomly sample 2.3M word pairs from the corpus of *Jeopardy!* clues available at j-archive.com (approximately 1% of the word pairs observed in the corpus). Using the GloVe Wikipedia 2014 and Gigaword 5 embeddings, I compute the cosine similarity, the Euclidean norm of the difference vector, and the log-rank frequency ratio in the GloVe corpus. The figure below plots the functional relationship between these quantities (linear $R^2 = 0.6369$)

```{r, echo=F, message=F, fig.height=4, fig.width=4.5, fig.align="center"}
setwd("/Users/akindel/Code/jeopardy")

library(here)

jtok_dists <- read_csv(here("data", "jtok_glove_cs_13July21.csv"))

jtok_dists %>%
  slice_sample(prop=0.01) %>%
  mutate(rnr = log(rn / rn2),
         rnr_b = log(rn * rn2)) %>%
  rowwise() %>%
  mutate(d_norm = norm(glove_d[rn,] - glove_d[rn2,], type="2")) ->
  test

test %>%
  ggplot(aes(x=cos_sim, y=d_norm, color=abs(rnr))) +
  geom_point()
```

Consider the case where $A, B$ and $Y$ are fixed. The difference vector $T$ is superimposed on the chord $T^*$ exactly when $X$ has magnitude $||X|| = ||Y||$. Perturbations away from $||X|| = ||Y||$ (or $||A|| = ||B||$) induce frequency bias when the underlying quantities are compared. This can be understood in terms of the difference between $T$ and $T^*$. When $||X|| > ||Y||$, its projection onto the unit circle, $\hat{X}$, will tend to be further from $Y$ than the point where the chord $T^*$ begins, $X^*$. The opposite holds when $||X|| < ||Y||$; the projection point $\hat{X}$ will tend to be closer to $Y$ than the chord endpoint $X^*$.
- TODO: Is this always true? Does it depend on the linear ordering of the orientations of the vectors relative to each other?
  
  When the ratio $||B||/||X||$ is very large relative to $||A||/||X||$, the comparison approaches the angle between the  diameter ending at $A$ and the tangent line at $Y$. In this regime, when $A \perp Y$, $T$ and $D$ are parallel, and the weight is 1; when $A$ and $Y$ are parallel, $T$ and $D$ are perpendicular, and the weight is 0. When these magnitude ratios are equalized, the range of potential weights is constrained away from these extremes.


The issue is not with the word embeddings -- which are merely more or less honest about the role of scale as a dimension of linguistic meaning -- but the attempt to apply compare cosine similarities with two different geometries using standard arithmetic, and the identification of this estimand with bias. The underlying word embedding model is in a sense innocent, not because vector space models do not encode scale (they do) but because the way the question has been phrased mathematically is not itself invariant to scale.

The use of pre-defined word lists is also problematic because it constrains our attention to a small and unrepresentative region of the word embedding model, and presumes that this is an adequate set of items to cover the concept we are interested in measuring. This plays on our sense that the words in each list obviously belong to the category in some general sense, leading us to overlook the relatively arbitrary mathematical structure of the lists: the original authors selected names that maximize their frequency ratio conditional on race (Greenwald et al. 1998; Bertrand \& Mullainathan 2004; Caliskan, Bryson & Narayanan 2017) without considering the contribution of the absolute frequencies of the names to this ratio.



An important but relatively uninteresting takeaway is that the WEAT measure is very sensitive to the relative lengths of the word lists. This is uninteresting in the sense that researchers can for the most part address this problem by using four word lists that are the same length. A slightly more interesting takeaway (which is, alas, uninteresting for the same reason) is that randomizing the order of the lists is important, because the original ordering of the lists can smuggle scale information into the paired comparisons WEAT operates on. But most interesting is that equalizing the lengths of the word lists and randomizing their order does not completely solve the problem, because the lists imply different relative term frequency distributions in different parts of the vector space.


```{r, echo=F, message=F, fig.height=4, fig.width=4, fig.align="center", fig.cap=""}
library(tidyverse)
library(latex2exp)
theme_set(theme_bw())
test <- data.frame(x=seq(-1, 1, by=0.005), y=seq(-1, 1, by=0.005)) %>% expand(x, y) %>% mutate(z=x + y)
test %>%
  ggplot(aes(x, y, fill=z)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(x=latex2exp::TeX("$cos(\\theta_{WA})$"),
       y=latex2exp::TeX("$cos(\\theta_{WB})$"),
       fill=latex2exp::TeX("$cos(\\theta_{WA}) + cos(\\theta_{WB})$")) +
  theme(legend.position='bottom')
```
||D^+_{abc}|| &= \sqrt{3 + 2\cos(a, b) + 2\cos(b, c) + 2\cos(a, c)}\\




An essential yet often overlooked aspect of cosine similarity is its connection to the Pearson product-moment correlation coefficient $r(a, b)$. In particular, like this closely related quantity, the bounded range of the cosine ratio implies that their arithmetic mean is structurally related to their distribution, resulting in a sampling distribution that is skewed by itself. It is well-established that the mean of these scaled distributions is asymptotically normal, and the extensive literature on correlation discusses normalizing transformations to overcome this issue (Fisher 1928; Hotelling 1953; Meng, Rosenthal & Rubin 1992). 


Jolicoeur, P. and J.E. Mosimann. 1960. "Size and shape variation in the painted turtle: A principal component analysis." *Growth* 24(4): 339-354.
Cvetković, D. M., and I. Gutman. 1977. "Note on branching." *Croatica Chemica Acta* 49(1): 115-121.

Carpenter, Jr., J. A., H. A. Fitzhugh, T. C. Cartwright, R. C. Thomas, and A. A. Melton. 1978. "Principal Components for Cow Size and Shape." *Journal of Animal Science* 46(2): 370-375.

Balaban, A.T., I. Motoc, D. Bonchev, and O. Mekenyan. 1983. "Topological indices for structure-activity correlations." In V. Austel, A. T. Balaban, D. Bonchev, M. Charton, T. Fujita, H. Iwamura, O. Mekenyan, and I. Motoc (eds.), *Steric Effects in Drug Design*: 21-55.





# # wgs <- t(apply(weat_white_glove, 1, scale2))[sample(1:nrow(weat_white_glove)),]  # Scaled 
# wgs_rand_r400000 <- t(apply(glove_d[sample(1:400000, 1000),], 1, scale2))  # Scaled 
# wgs_rand_r100000 <- t(apply(glove_d[sample(1:100000, 1000),], 1, scale2))  # Scaled 
# wgs_rand_r20000 <- t(apply(glove_d[sample(1:20000, 1000),], 1, scale2))  # Scaled 
# wgs_rand_r5000 <- t(apply(glove_d[sample(1:5000, 1000),], 1, scale2))  # Scaled 
# 
# # Plot a single term against an expanding sum vector of random terms from embedding (rarity stratified)
# # Note varying X and Y axes.
# test_term <- "silver"
# expand.grid(i=2:500) %>%
#   rowwise() %>%
#   mutate(dnorm_r5000 = norm(colSums(wgs_rand_r5000[1:i,]), "2"),
#          dnorm_r20000 = norm(colSums(wgs_rand_r20000[1:i,]), "2"),
#          dnorm_r100000 = norm(colSums(wgs_rand_r100000[1:i,]), "2"),
#          dnorm_r400000 = norm(colSums(wgs_rand_r400000[1:i,]), "2"),
#          sim_r5000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r5000[1:i,])),
#          sim_r20000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r20000[1:i,])),
#          sim_r100000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r100000[1:i,])),
#          sim_r400000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r400000[1:i,]))) %>%
#   mutate(cmean_sc = mean(.$sim_r5000[1:i-1]),
#          cmean_sm = mean(.$sim_r20000[1:i-1]),
#          cmean_sr = mean(.$sim_r100000[1:i-1]),
#          cmean_sv = mean(.$sim_r400000[1:i-1])) ->
#   term_simd
# 
# # Plot CS(w_i, D) against ||D|| for arbitrary choice of w_i
# # The behavior of the function is really ugly
# # You can plot the mean by dividing dnorm_x by i (the number of observations in the sum vector)
# # The behavior of the mean is really ugly compared to the sum.
# grid.arrange(term_simd %>%
#                ggplot(aes(x=dnorm_r5000, y=sim_r5000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r5000, y=cmean_sc, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r5000)), linetype="dashed") +
#                labs(title=latex2exp::TeX("proj($v_{Alex}$, $D_{rand})$")) +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r20000, y=sim_r20000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r20000, y=cmean_sm, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r20000)), linetype="dashed") +
#                labs(title=latex2exp::TeX("proj($v_{Alex}$, $D_{rand})$")) +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r100000, y=sim_r100000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r100000, y=cmean_sr, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r100000)), linetype="dashed") +
#                labs(title=latex2exp::TeX("proj($v_{Alex}$, $D_{rand})$")) +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r400000, y=sim_r400000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r400000, y=cmean_sv, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r400000)), linetype="dashed") +
#                labs(title=latex2exp::TeX("proj($v_{Alex}$, $D_{rand})$")) +
#                theme(legend.position="bottom"),
#              ncol=2)



# As a function of the major and minor diagonals separately
# Also interesting to look at just ab_cs or just anorm*bnorm
# random_angles %>% ggplot(aes(y=anorm*bnorm*ab_cs, x=plus_dnorm)) + geom_point()
# random_angles %>% ggplot(aes(y=anorm*bnorm*ab_cs, x=minus_dnorm)) + geom_point()

# Plot the inner product
# You can see how the angular decomposition classifies the parallelograms
# random_angles %>%
#   ggplot(aes(x=anorm*bnorm, y=ab_cs)) +
#   geom_point(aes(color=plus_dnorm, size=minus_dnorm)) +
#   geom_smooth(color="grey75") +
#   scale_color_viridis_c()

# An interesting thing to do is to pull out the comparisons by norm product weight
# There's a pretty clear linear relationship in the norm product weight in the close ones but not the far ones
# epsil <- 0.1
# random_angles_50 %>%
#     filter(abs(anorm - bnorm) < epsil) %>%
#     ggplot(aes(x=anorm*bnorm, y=ab_cs)) +
#     geom_point(aes(color=plus_dnorm, size=minus_dnorm)) +
#     geom_smooth(color="grey75") +
#     scale_color_viridis_c()




# Plot CS(w_i, D) against ||D|| for arbitrary choice of w_i
# The behavior of the function is really ugly
# You can plot the mean by dividing dnorm_x by i (the number of observations in the sum vector)
# The behavior of the mean is really ugly compared to the sum.
# grid.arrange(term_simd %>%
#                ggplot(aes(x=dnorm_r5000, y=sim_r5000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r5000, y=cmean_sc, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r5000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 5k terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r20000, y=sim_r20000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r20000, y=cmean_sm, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r20000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 5% terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r40000, y=sim_r40000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r40000, y=cmean_sk, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r40000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 10% terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r80000, y=sim_r80000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r80000, y=cmean_sl, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r80000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 20% terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r100000, y=sim_r100000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r100000, y=cmean_sr, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r100000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 25% terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r200000, y=sim_r200000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r200000, y=cmean_sp, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r200000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (top 50% terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              term_simd %>%
#                ggplot(aes(x=dnorm_r400000, y=sim_r400000, color=i)) +
#                geom_point() +
#                geom_line(aes(x=dnorm_r400000, y=cmean_sv, color=i)) +
#                geom_hline(aes(yintercept=mean(sim_r400000)), linetype="dashed") +
#                labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
#                     x="||D|| (all terms)", y="cos(v, D)") +
#                theme(legend.position="bottom"),
#              ncol=2)




# Plot a single term against an expanding sum vector of random terms from embedding (rarity stratified)
# Note varying X and Y axes.
k <- 100
expand.grid(i=2:k) %>%
  rowwise() %>%
  mutate(# dnorm_r5000 = norm(colSums(wgs_rand_r5000[1:i,]), "2"),
    # dnorm_r20000 = norm(colSums(wgs_rand_r20000[1:i,]), "2"),
    dnorm_r40000 = norm(colSums(wgs_rand_r40000[1:i,]), "2"),
    # dnorm_r80000 = norm(colSums(wgs_rand_r80000[1:i,]), "2"),
    dnorm_r100000 = norm(colSums(wgs_rand_r100000[1:i,]), "2"),
    dnorm_r200000 = norm(colSums(wgs_rand_r200000[1:i,]), "2"),
    dnorm_r400000 = norm(colSums(wgs_rand_r400000[1:i,]), "2"),
    # sim_r5000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r5000[1:i,])),
    # sim_r20000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r20000[1:i,])),
    sim_r40000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r40000[1:i,])),
    # sim_r80000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r80000[1:i,])),
    sim_r100000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r100000[1:i,])),
    sim_r200000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r200000[1:i,])),
    sim_r400000 = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_r400000[1:i,])),
    # simsum_r5000 = sum_cossim(i, wgs_rand_r5000),
    # simsum_r20000 = sum_cossim(i, wgs_rand_r20000),
    simsum_r40000 = sum_cossim(i, wgs_rand_r40000),
    # simsum_r80000 = sum_cossim(i, wgs_rand_r80000),
    simsum_r100000 = sum_cossim(i, wgs_rand_r100000),
    simsum_r200000 = sum_cossim(i, wgs_rand_r200000),
    simsum_r400000 = sum_cossim(i, wgs_rand_r400000),
    simmean_r40000 = simsum_r40000/i,
    simmean_r100000 = simsum_r100000/i,
    simmean_r200000 = simsum_r200000/i,
    simmean_r400000 = simsum_r400000/i) %>%
  mutate(# cmean_sc = mean(.$sim_r5000[1:i-1]),
    # cmean_sm = mean(.$sim_r20000[1:i-1]),
    cmean_sk = mean(.$sim_r40000[1:i-1]),
    # cmean_sl = mean(.$sim_r80000[1:i-1]),
    cmean_sr = mean(.$sim_r100000[1:i-1]),
    cmean_sp = mean(.$sim_r200000[1:i-1]),
    cmean_sv = mean(.$sim_r400000[1:i-1])) ->
  term_simd

grid.arrange(term_simd %>%
               # ggplot(aes(x=simsum_r200000, y=sim_r200000, color=i)) +
               ggplot(aes(color=dnorm_r40000, y=simsum_r40000, x=i)) +
               geom_point() +
               # geom_line(aes(x=simsum_r40000, y=cmean_sk, color=i)) +
               # geom_line(aes(x=dnorm_r40000, y=cmean_sk, color=i)) +
               # geom_hline(aes(yintercept=mean(sim_r40000)), linetype="dashed") +
               labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
                    color="||D|| (top 10% terms)", y=latex2exp::TeX("$\\sum_{d \\in D} cos(v, d)$")) +
               theme(legend.position="bottom"),
             term_simd %>%
               # ggplot(aes(x=simsum_r200000, y=sim_r200000, color=i)) +
               ggplot(aes(color=dnorm_r100000, y=simsum_r100000, x=i)) +
               geom_point() +
               # geom_line(aes(x=simsum_r100000, y=cmean_sr, color=i)) +
               # geom_line(aes(x=dnorm_r100000, y=cmean_sr, color=i)) +
               # geom_hline(aes(yintercept=mean(sim_r100000)), linetype="dashed") +
               labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
                    color="||D|| (top 25% terms)", y=latex2exp::TeX("$\\sum_{d \\in D} cos(v, d)$")) +
               theme(legend.position="bottom"),
             term_simd %>%
               # ggplot(aes(x=simsum_r200000, y=sim_r200000, color=i)) +
               ggplot(aes(color=dnorm_r200000, y=simsum_r200000, x=i)) +
               geom_point() +
               # geom_line(aes(x=simsum_r200000, y=cmean_sp, color=i)) +
               # geom_line(aes(x=dnorm_r200000, y=cmean_sp, color=i)) +
               # geom_hline(aes(yintercept=mean(sim_r200000)), linetype="dashed") +
               labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
                    color="||D|| (top 50% terms)", y=latex2exp::TeX("$\\sum_{d \\in D} cos(v, d)$")) +
               theme(legend.position="bottom"),
             term_simd %>%
               # ggplot(aes(x=simsum_r200000, y=sim_r200000, color=i)) +
               ggplot(aes(color=dnorm_r400000, y=simsum_r400000, x=i)) +
               geom_point() +
               # geom_line(aes(x=simsum_r400000, y=cmean_sv, color=i)) +
               # geom_line(aes(x=dnorm_r400000, y=cmean_sv, color=i)) +
               # geom_hline(aes(yintercept=mean(sim_r400000)), linetype="dashed") +
               labs(title=latex2exp::TeX(paste0("proj($v_{", test_term, "}$, $D_{rand})$")),
                    color="||D|| (all terms)", y=latex2exp::TeX("$\\sum_{d \\in D} cos(v, d)$")) +
               theme(legend.position="bottom"),
             ncol=2)
}

Tversky, A. 1977. "Features of similarity." *Psychological Review* 84(4): 327-352.


# The key feature of this quantity is its reduction of the dimensionality of the comparison between the vectors. While $A$ and $B$ each have $p$ entries, they are also always coplanar, so the angle $\theta_{AB}$ formed between $A$ and $B$ at the origin along this plane is a scalar. This feature is responsible for the ratio's reputation as a scale-independent measure of relatedness: 

This property of the  that make it appealing for applied research, particularly in high-dimensional data analysis settings. In particular, it is prized as a measure of similarity that is scale-independent; that is, it does not depend on the respective magnitudes of the two vectors under comparison. This property makes it fast to compute

Connor and Ross (2012) observed that there has been relatively little direct investigation of the statistical properties of cosine similarity, particularly as the number of underlying dimensions in the vector space increases. But the cosine similarity is equivalent to the Pearson product-moment correlation coefficient, and the necessary results have been established on more solid algebraic footing in the statistical literature on correlation estimation (Fisher 1928; Hotelling 1953). The sampling distribution of the sample correlation coefficient is non-normal; it is skewed in a way that depends on the shape of the dataset ().

Hotelling, H. 1953. "New light on the correlation coefficient and its transforms." *Journal of the Royal Statistical Society. Series B (Methodological)* 15(2): 193-232.
Fisher, R. A. 1928. "The general sampling distribution of the multiple correlation coefficient." *Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character* 121(788), 654-673.
Bai, Z. and J.W. Silverstein. 2010. *Spectral Analysis of Large Dimensional Random Matrices.* New York: Springer.

# Let y be a stochastic function of two random bivariate vectors w. normal components having slightly different means
# Two common such functions: the sum vector [x1+x2, x3+x4] and the mean vector [(x1+x3)/2, (x2+x4)/2]
# And then consider two smooth additive regressions, y ~ x1+x2+x3+x4 and y ~ x1a+x2a+x3a+x4a
# 
i <- 5000
data.frame(x1=rnorm(i, mean=0.4, sd=2),
           x2=rnorm(i, mean=-0.6, sd=2),
           x3=rnorm(i, mean=0.1, sd=2),
           x4=rnorm(i, mean=0.3, sd=2)) %>%
  mutate(y1=x1+x2+x3+x4+rnorm(i, mean=1, sd=2),
         y2=(x1+x3)/2 + (x2+x4)/2 + ((x1+x3)/2 * (x2+x4)/2) + rnorm(i, mean=1, sd=2))  %>%
  rowwise() %>%
  mutate(n12 = norm(c(x1, x2), "2"),
         n34 = norm(c(x3, x4), "2"),
         abcs = lsa::cosine(c(x1, x2), c(x3, x4)),
         x1a=x1/n12,
         x2a=x2/n12,
         x3a=x3/n34,
         x4a=x4/n34) ->
  dft


# Do a bunch of paired R3+ space comparisons and look at their common embedding space metrics
# They're all totally coincident in the row space but vary smoothly in the column space
k <- 3
rmax <- 40000
glovwind <- glove_tokens[1:rmax]
data.frame(i=1:10) %>%
  rowwise() %>%
  mutate(ktr = krz_trace(sample(glovwind, k), sample(glovwind, k))) ->
  test




# WEAT name range sample
# Look at 100 random tokens in the WEAT name frequency ranges
# wbt <- which(etokens %in% weat_black)
# wwt <- which(etokens %in% weat_white)
# random_angles_100_uneven <- make_angles_2(100, a_rmin=min(wwt), a_rmax=max(wwt), b_rmin=min(wbt), b_rmax=max(wbt))
# p <- plot_angle_manifold(random_angles_100_uneven)


# Function to plot the inner product space instead
# 
# plot_angle_manifold_ip <- function(r_angles) {
#   r_angles %>%
#     ggplot(aes(x=anorm*bnorm, y=ab_cs)) +
#     geom_point(aes(color=ab_ip, size=ab_ip)) +
#     # geom_rug(aes(color=plus_dnorm), outside=T) +
#     geom_smooth(method="loess", color="grey75") +
#     scale_color_viridis_c() +
#     theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
#     labs(x=latex2exp::TeX("$||A||*||B||$"), y=latex2exp::TeX("$cos(A, B)$"),
#          color=latex2exp::TeX("||D(+)||"), size=latex2exp::TeX("||D(-)||"))
# }








# Extra: Cosine as a mixture of distributions of related and unrelated terms

# This is almost always suggestive of a Gaussian mixture (related and unrelated terms)
# This can be modeled explicitly
t1t_2g <- mixtools::normalmixEM(test_1term$angle)

# Plot histogram
# Superimpose mixture of Gaussians, optionally
mean_angle <- mean(test_1term$angle, na.rm=T)
median_angle <- median(test_1term$angle, na.rm=T)
bw <- 1
test_1term %>%
  ggplot(aes(x=angle)) +
  geom_histogram(aes(y=..density..), binwidth=bw) +
  geom_vline(xintercept=mean_angle, linetype="dashed", color="dodgerblue") +
  geom_vline(xintercept=median_angle, linetype="dashed", color="tomato") +
  geom_vline(xintercept=90, linetype="dotted") +
  # stat_function(fun=dnorm, args=list(mean=t1t_2g$mu, sd=t1t_2g$sigma), n=1000, alpha=0.3) +
  scale_x_continuous(limits=c(0, 120), breaks=seq(0, 120, by=10))


# We are often looking at sum vectors or mean vectors in comparison to w_j
# What is the distribution of angles between the implied subspaces?
augment_tri <- function(embm, ratri) {
  ratri %>%
    rowwise() %>%
    mutate(bcnorm = norm(embm[bterm,] + embm[cterm,], "2"),
           a.bc_cs = lsa::cosine(as.numeric(embm[aterm,]), as.numeric(embm[bterm,] + embm[cterm,]))[1],
           a.bc_ip = as.numeric(embm[aterm,]) %*% as.numeric(embm[bterm,] + embm[cterm,]),
           bc.basis.ang = lsa::cosine(svd(embm[c(bterm, cterm),])$v[,1],
                                      svd(rbind(embm[bterm,]/bnorm, embm[cterm,]/cnorm))$v[,1]),
           a.bc.basis.ang = lsa::cosine(svd(rbind(embm[aterm,], embm[bterm,]+embm[cterm,]))$v[,1],
                                        svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm + embm[cterm,]/cnorm))$v[,1])) ->
    ratri.aug
  return(ratri.aug)
}

gt.nemb.aug <- augment_tri(nemb, gt.nemb)

# Compare the normalization column deflections of BC-BC' and ABC-ABC'
# When the BC distortion is larger, the ABC normalization is less well-aligned
gt.nemb.aug %>%
  ggplot(aes(x=abs(bc.basis.ang), y=abs(a.bc.basis.ang), color=a.bc_cs)) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom")

gt.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(ab.pr.ang * ac.pr.ang, 12)) %>%
  arrange(ab.pr.ang * ac.pr.ang) %>%
  ggplot(aes(x=abs(abc.ang), y=abs(abc.ang.sc), color=ab.pr.ang * ac.pr.ang)) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile)



# Construct full LNPW and root log frequency product distributions, ~856M comparisons
# This can take a long time since we're looking at the better part of 1B term pairs.

tfn %>%
  ungroup() %>%
  dplyr::select(t1=feature, t2=feature) %>%
  tidyr::expand(t1, t2) %>%
  filter(t1 != t2) %>%
  left_join(tfn %>% dplyr::select(t1=feature, t1n=snorm, t1f=frequency), by="t1") %>%
  left_join(tfn %>% dplyr::select(t2=feature, t2n=snorm, t2f=frequency), by="t2") %>%
  mutate(lnpw=t1n*t2n, slfr=sqrt(log(t1f))*sqrt(log(t2f))) ->
  freqbiasmat




# Arora et al. plot, but for the local norm product weight; i.e. this is the same plot
#  but for the pair-directed inner product and not the self-directed inner product.
lnpw_plot <- function(fbm, subsample=1) {
  fbm %>%
    sample_frac(size=subsample) %>%
    ggplot(aes(x=sqrt(log(a$frequency))*sqrt(log(b$frequency)), y=nprod)) +
    geom_point() +
    geom_smooth()
}




\epigraph{\itshape Language is a labyrinth of paths. You approach from one side and know your way about; you approach the same place from another side and no longer know your way about.}{Ludwig Wittgenstein, \textit{Philosophical Investigations} \S 203 (1953)}


# linear_frequency_bias(make_dim_df(cs.male_black_polity, cs.female_black_polity, "male_black_polity", "female_black_polity", mv_male_black, mv_female_black, polity_words), "Black Male x Black Female (polity)")
# linear_frequency_bias(make_dim_df(cs.male_black_economy, cs.female_black_economy, "male_black_economy", "female_black_economy", mv_male_black, mv_female_black, economy_words), "Black Male x Black Female (economy)")
# linear_frequency_bias(make_dim_df(cs.male_black_culture, cs.female_black_culture, "male_black_culture", "female_black_culture", mv_male_black, mv_female_black, culture_words), "Black Male x Black Female (culture)")
# linear_frequency_bias(make_dim_df(cs.male_black_domestic, cs.female_black_domestic, "male_black_domestic", "female_black_domestic", mv_male_black, mv_female_black, domestic_words), "Black Male x Black Female (domestic)")
# 
# # White male <-> White female
# linear_frequency_bias(make_dim_df(cs.male_white_polity, cs.female_white_polity, "male_white_polity", "female_white_polity", mv_male_white, mv_female_white, polity_words), "White Male x White Female (polity)")
# linear_frequency_bias(make_dim_df(cs.male_white_economy, cs.female_white_economy, "male_white_economy", "female_white_economy", mv_male_white, mv_female_white, economy_words), "White Male x White Female (economy)")
# linear_frequency_bias(make_dim_df(cs.male_white_culture, cs.female_white_culture, "male_white_culture", "female_white_culture", mv_male_white, mv_female_white, culture_words), "White Male x White Female (culture)")
# linear_frequency_bias(make_dim_df(cs.male_white_domestic, cs.female_white_domestic, "male_white_domestic", "female_white_domestic", mv_male_white, mv_female_white, domestic_words), "White Male x White Female (domestic)")
# 
# # White female <-> Black female
# linear_frequency_bias(make_dim_df(cs.female_white_polity, cs.female_black_polity, "female_white_polity", "female_black_polity", mv_female_white, mv_female_black, polity_words), "White Female x Black Female (polity)")
# linear_frequency_bias(make_dim_df(cs.female_white_economy, cs.female_black_economy, "female_white_economy", "female_black_economy", mv_female_white, mv_female_black, economy_words), "White Female x Black Female (economy)")
# linear_frequency_bias(make_dim_df(cs.female_white_culture, cs.female_black_culture, "female_white_culture", "female_black_culture", mv_female_white, mv_female_black, culture_words), "White Female x Black Female (culture)")
# linear_frequency_bias(make_dim_df(cs.female_white_domestic, cs.female_black_domestic, "female_white_domestic", "female_black_domestic", mv_female_white, mv_female_black, domestic_words), "White Female x Black Female (domestic)")
# 
# # Black male <-> White male
# linear_frequency_bias(make_dim_df(cs.male_black_polity, cs.male_white_polity, "male_black_polity", "male_white_polity", mv_male_black, mv_male_white, polity_words), "Black Male x White Male (polity)")
# linear_frequency_bias(make_dim_df(cs.male_black_economy, cs.male_white_economy, "male_black_economy", "male_white_economy", mv_male_black, mv_male_white, economy_words), "Black Male x White Male (economy)")
# linear_frequency_bias(make_dim_df(cs.male_black_culture, cs.male_white_culture, "male_black_culture", "male_white_culture", mv_male_black, mv_male_white, culture_words), "Black Male x White Male (culture)")
# linear_frequency_bias(make_dim_df(cs.male_black_domestic, cs.male_white_domestic, "male_black_domestic", "male_white_domestic", mv_male_black, mv_male_white, domestic_words), "Black Male x White Male (domestic)")
# 
# # White male <-> Black female
# linear_frequency_bias(make_dim_df(cs.male_white_polity, cs.female_black_polity, "male_white_polity", "female_black_polity", mv_male_white, mv_female_black, polity_words), "White Male x Black Female (polity)")
# linear_frequency_bias(make_dim_df(cs.male_white_economy, cs.female_black_economy, "male_white_economy", "female_black_economy", mv_male_white, mv_female_black, economy_words), "White Male x Black Female (economy)")
# linear_frequency_bias(make_dim_df(cs.male_white_culture, cs.female_black_culture, "male_white_culture", "female_black_culture", mv_male_white, mv_female_black, culture_words), "White Male x Black Female (culture)")
# linear_frequency_bias(make_dim_df(cs.male_white_domestic, cs.female_black_domestic, "male_white_domestic", "female_black_domestic", mv_male_white, mv_female_black, domestic_words), "White Male x Black Female (domestic)")
# 
# # White female <-> Black male
# linear_frequency_bias(make_dim_df(cs.female_white_polity, cs.male_black_polity, "female_white_polity", "male_black_polity", mv_female_white, mv_male_black, polity_words), "White Female x Black Male (polity)")
# linear_frequency_bias(make_dim_df(cs.female_white_economy, cs.male_black_economy, "female_white_economy", "male_black_economy", mv_female_white, mv_male_black, economy_words), "White Female x Black Male (economy)")
# linear_frequency_bias(make_dim_df(cs.female_white_culture, cs.male_black_culture, "female_white_culture", "male_black_culture", mv_female_white, mv_male_black, culture_words), "White Female x Black Male (culture)")
# linear_frequency_bias(make_dim_df(cs.female_white_domestic, cs.male_black_domestic, "female_white_domestic", "male_black_domestic", mv_female_white, mv_male_black, domestic_words), "White Female x Black Male (domestic)")
# 
# # Authority (Fig. 5)
# linear_frequency_bias(make_dim_df(cs.male_black_authority, cs.female_black_authority, "male_black_authority", "female_black_authority", mv_male_black, mv_female_black, authority_words), "Black Male x Black Female (authority)")
# linear_frequency_bias(make_dim_df(cs.male_white_authority, cs.female_white_authority, "male_white_authority", "female_white_authority", mv_male_white, mv_female_white, authority_words), "White Male x White Female (authority)")
# linear_frequency_bias(make_dim_df(cs.female_white_authority, cs.female_black_authority, "female_white_authority", "female_black_authority", mv_female_white, mv_female_black, authority_words), "White Female x Black Female (authority)")
# linear_frequency_bias(make_dim_df(cs.male_black_authority, cs.male_white_authority, "male_black_authority", "male_white_authority", mv_male_black, mv_male_white, authority_words), "Black Male x White Male (authority)")
# linear_frequency_bias(make_dim_df(cs.male_white_authority, cs.female_black_authority, "male_white_authority", "female_black_authority", mv_male_white, mv_female_black, authority_words), "White Male x Black Female (authority)")
# linear_frequency_bias(make_dim_df(cs.female_white_authority, cs.male_black_authority, "female_white_authority", "male_black_authority", mv_female_white, mv_male_black, authority_words), "White Female x Black Male (authority)")



The approach I take in this paper is modeled after the literature on regression diagnostics and residual analysis (Cook & Weisberg 1982; Belsley, Kuh & Welsch 1980). Cook and Weisberg write that diagnostic statistics serve two purposes: “they may result in the recognition of important phenomena that might otherwise have gone unnoticed” and they “can be used to suggest appropriate remedial action to the analysis of the model“ (p. 2). I focus on the former task; the latter task I leave for future work, although I comment on some potential directions at the end of the paper. My aim is to enable researchers in the social sciences to notice something about word embedding models that has been “hidden in plain sight” (Zerubavel 2015), and I hope that doing so will provoke further work on ways of measuring aggregate similarity over embedding vector subspaces.


plot_fbias <- function(dimdf, fbtitle) {
  dimdf %>%
    ggplot(aes(x=nprod, y=cs, color=identity)) +
    geom_point() +
    geom_hline(yintercept=mean(dimdf %>% filter(identity == first(identity)) %$% cs),
               color="#00BFC4", linetype="dashed") +
    geom_hline(yintercept=mean(dimdf %>% filter(identity == last(identity)) %$% cs),
               color="#F8766D", linetype="dashed") +
    geom_smooth(method="gam") +
    ggtitle(fbtitle)
}



## Dependence in pair-structured cosine similarity matrices

```{r arithmetic_strats, echo=F, out.width="\\textwidth", fig.cap="Arithmetic comparisons on $G(A, B)$; matrix entries represent $\\cos(A_i, B_j)$."}
knitr::include_graphics("./arithmetic_comparisons.png")
```

Comparisons (1) and (3) involve reusing vectors in multiple comparisons. Vectors are reused in multiple comparisons --> bilinear clustering of cosine similarities. 

The cosine similarity between two word vectors is exactly equivalent to the Pearson product-moment correlation coefficient between them (Fisher XXXX; Hotelling 1953; Dunn \& Clark 1971; Rodgers \& Nicewander 1980; Steiger 1980). The distribution of this quantity is known to be skewed at small n, large |rho|. In general there are further problems with aggregation the larger the cosine matrix becomes. See Steiger (1980) on summarizing correlation matrices with paired dependence structure. This literature has implications for double-dipping on terms to compute multiple angles to the rest of the space from a single position, and the need to account for this when estimating uncertainty.

The arithmetic mean sample cosine similarity is generally non-zero in the full embedding space and in random term subspaces, and it is (in part) a function of the word frequency distribution of the component vectors. Additionally, it is unclear whether there is a well-defined notion of population cosine similarity in natural language (Baayen 2001, ch. 10; Piantadosi 2014).


# 
# polity_svd <- svd(polity_vmat)
# economy_svd <- svd(economy_vmat)
# culture_svd <- svd(culture_vmat)
# domestic_svd <- svd(domestic_vmat)
# authority_svd <- svd(authority_vmat)
# 
# # # Normalized versions
# wordspace::normalize.rows(polity_vmat) -> polity_vn
# wordspace::normalize.rows(economy_vmat) -> economy_vn
# wordspace::normalize.rows(culture_vmat) -> culture_vn
# wordspace::normalize.rows(domestic_vmat) -> domestic_vn
# wordspace::normalize.rows(authority_vmat) -> authority_vn
# 
# 
# # Additive race/gender subspaces
# wbf <- expand_grid(fwl_black, fwl_female)
# wwf <- expand_grid(fwl_white, fwl_female)
# wbm <- expand_grid(fwl_black, fwl_male)
# wwm <- expand_grid(fwl_white, fwl_male)
# bf <- as.matrix(nemb[wbf$fwl_black,] + nemb[wbf$fwl_female,])
# wf <- as.matrix(nemb[wwf$fwl_white,] + nemb[wwf$fwl_female,])
# bm <- as.matrix(nemb[wbm$fwl_black,] + nemb[wbm$fwl_male,])
# wm <- as.matrix(nemb[wwm$fwl_white,] + nemb[wwm$fwl_male,])
# 
# bfs <- svd(bf)
# wfs <- svd(wf)
# bms <- svd(bm)
# wms <- svd(wm)


# Crone and Crosby (1995 Technometrics) provide a distance metric between column subspaces
#  based on Krzanowski 1979: sqrt(p - tr(v_Av_B)) where v_M is the column span of M.
ccdist <- function(A, B) {
  p <- ncol(A)
  Ap <- A %*% solve(t(A)%*%A) %*% t(A)
  Bp <- B %*% solve(t(B)%*%B) %*% t(B)
  prtr <- sum(diag(Ap %*% Bp))
  d <- sqrt(p - prtr)
  # return(d)
  return(prtr)
}

# ccdist(polity_svd$u, economy_svd$u)
# ccdist(wms$v, bfs$v)

# Anisotropy interlude
# 
# # An interesting plot: the variance-covariance matrix of these subspaces
# # data.frame(v1 = colnames(economy_vmat), data.frame(cov(economy_vmat))) %>%
# data.frame(v1 = colnames(nemb), data.frame(cov(nemb))) %>%
#   pivot_longer(starts_with("X")) %>%
#   rename(v2=name) %>%
#   ggplot(aes(x=v1, y=v2, fill=value)) +
#   geom_tile() +
#   scale_fill_viridis_c()
# 
# # An important problem is that this methodology selects subspaces with different
# #  covariance structures from each other and the global vector space. This can be
# #  demonstrated by showing that the subspaces have different spectra.
# 
# # Mahalanobis distance between sets.
# maha_econ_dom <- data.frame(economy_words, md=mahalanobis(economy_vmat, colMeans(domestic_vmat), pracma::pinv(cov(domestic_vmat)), inverted=T))
# maha_dom_econ <- data.frame(domestic_words, md=mahalanobis(domestic_vmat, colMeans(economy_vmat), pracma::pinv(cov(economy_vmat)), inverted=T))
# maha_pol_econ <- data.frame(polity_words, md=mahalanobis(polity_vmat, colMeans(economy_vmat), pracma::pinv(cov(economy_vmat)), inverted=T))
# 
# # Hausdorff distance from subspace A to subspace B is the maximum distance
# #  from a vector in subspace A to its closest vector in subspace B.
# #  The quantity is non-commutative when the subspaces aren't the same shape;
# #   this is another way of saying the distance matrix is asymmetric.
# #  You can use the difference as a measure of anisotropy
# #   - you can perturb this to look at how sensitive the anisotropy is to leave-one-out
# polity_economy_dm <- as.matrix(pdist::pdist(polity_vmat, economy_vmat))
# #polity_economy_dm %>% reshape2::melt() %>% ggplot(aes(x=Var1, y=Var2, fill=value)) + geom_tile()  # Asymmetric
# hd_polity_economy_left <- max(apply(polity_economy_dm, 1, min))  # DH(polity->economy)
# hd_polity_economy_right <- max(apply(polity_economy_dm, 2, min))  # DH(economy->polity)
# 
# # There is often a small set of target points that induces the metric for a large
# #  set of the starting points, in both directions. The modal closest point tends to be
# #  the maximally far point in the set of closest points (i.e. the Hausdorff inducing point).
# hd_polity_economy_left.v <- apply(polity_economy_dm, 1, which.min)
# hd_polity_economy_right.v <- apply(polity_economy_dm, 2, which.min)
# hd_polity_economy_left.i <- as.integer(names(which.max(table(apply(polity_economy_dm, 1, which.min)))))
# hd_polity_economy_right.i <- as.integer(names(which.max(table(apply(polity_economy_dm, 2, which.min)))))
# hdpe_left.it <- economy_words[hd_polity_economy_left.i]
# hdpe_right.it <- polity_words[hd_polity_economy_right.i]
# 
# # Another example of anisotropic comparison
# # This one is less asymmetric
# culture_domestic_dm <- as.matrix(pdist::pdist(culture_vmat, domestic_vmat))
# hd_culture_domestic_left <- max(apply(culture_domestic_dm, 1, min))  # DH(culture->domestic)
# hd_culture_domestic_right <- max(apply(culture_domestic_dm, 2, min))  # DH(domestic->culture)
# hd_culture_domestic_left.v <- apply(culture_domestic_dm, 1, which.min)
# hd_culture_domestic_right.v <- apply(culture_domestic_dm, 2, which.min)
# 
# culture_economy_dm <- as.matrix(pdist::pdist(culture_vmat, economy_vmat))
# hd_culture_economy_left <- max(apply(culture_economy_dm, 1, min))  # DH(culture->economy)
# hd_culture_economy_right <- max(apply(culture_economy_dm, 2, min))  # DH(economy->culture)
# hd_culture_economy_left.v <- apply(culture_economy_dm, 1, which.min)
# hd_culture_economy_right.v <- apply(culture_economy_dm, 2, which.min)
# 


# Back to frequency bias...!


A weight-based interpretation of cosine similarity clarifies the role of the ratio in producing frequency-biased quantities from word embedding models. Taking the cosine of two word vectors enables us to estimate how far apart they are adjusted by their respective positions in the underlying word frequency distribution. However, the way that cosine similarity accomplishes this transformation is not constant everywhere in the original inner product space. The local normalization weight function encodes a notion of similarity that prefers rare words: holding the inner product constant, cosine similarity increases as the product of the vector norms becomes small. Consequently, conventional arithmetic summaries of cosine similarities tend to show strong correlations with measures of word frequency, particularly when the set of cosine similarities is determined by a pairwise correlation structure, as is common in applications of word embeddings in the social sciences. The root of frequency bias is an implicit assumption that the local *within-pair* interpretation of cosine similarity extends to the global *between-pair* mathematical characteristics of arithmetic summaries of cosine similarities. In practice, the amount of distortion implied by cosine similarity varies from application to application; I turn next to discussing a few cases of frequency bias in applied settings before returning to  



## Correlation structure and scale clustering

Another common strategy in the literature involves summarizing similarities between two sets of word vectors $X, Y$ of the form $\Psi_{X,Y}=\left\{\cos(X_i,Y_j)\right\}_{i=1, j=1}^{N,M}$. Each $X_i$ is reused in $j$ comparisons and each $Y_j$ is reused in $i$ comparisons. A closely related comparison within one set of vectors $X$ is $\Psi_{X,X}=\left\{\cos(X_i,X_j)\right\}_{i\neq j}^{N}$, where each vector is compared to the $N-1$ other vectors in the set.


In particular, although the exact shape of semantic vector spaces is somewhat model dependent, in general log-bilinear word embedding models are designed to estimate frequency-varying vector norms in order to optimally encode the underlying word-word co-frequency data (Mnih & Hinton 2008; Andreas & Klein 2015; Arora et al. 2016, 2018). Consequently, there is no "neutral" choice of weights that results in a distance metric that would allow us to ignore word frequency.^[A uniform distribution of meaning over frequency is not a neutral choice of weights; arguably, it also does not result in a realistic model, because we should avoid the assumption that the level of generality/specificity in discourse is constant.] Similarities are heuristic projections of semantic vector spaces that explicitly express the researcher's preference for associations with common or rare words in terms of a second normalized inner product space. They provide an *intentionally* distorted perspective on the inner product space that is useful for some purposes and less for others. Another implication of viewing similarity as a weight estimation procedure is that researchers need not apply any one similarity coefficient in order to construct informative models of variation in meaning using word embeddings. Many potential weighting schemes might be appropriate for different types of question, to the extent that they help researchers explore relevant regions of the original semantic vector space.

# Appendix B: Replications in progress


## “Aligning Differences: Discursive Diversity and Team Performance.” (Lix et al. 2022) 

\textit{In progress.}

## “The Automatic Analysis of Emotion in Political Speech Based on Transcripts.” (Cochrane et al. 2022)

\textit{In progress.}

The previous section suggests that cosine similarity can be understood as a non-Euclidean projection of the inner product space that distorts it in a way described by the local normalization weight function. This suggests the following geometric interpretation of the projection: it is smoothly bending the rare region of the underlying word frequency distribution toward the global mean vector. We can visualize this by looking at the mean cosine similarity at every frequency with the global mean vector (Figure \ref{fig:global_mean_csbend}). This brings our attention to the functional form of the cosine similarity operation. What enables this seemingly very simple quantity to perform such a complicated transformation of the inner product space?

\left<\Omega^{A}, X_i\right>&= \gamma_0\text{LNW}(\Omega^{A}, X_i) + \gamma_1\cos(\Omega^{B}, X_i)\text{LNW}(\Omega^{A}, X_i) + v_i \\
&+ \color{red}{\gamma_2 + \gamma_3\cos(\Omega^{B}, X_i) + \gamma_4\text{LNW}^{-1}(\Omega^{B}, X_i) + \gamma_5\left<\Omega^{B}, X_i\right>} \\
&+ \color{red}{\gamma_6\left<\Omega^{B}, X_i\right> \text{LNW}(\Omega^{A}, Y_i) + \gamma_7\frac{\text{LNW}(\Omega^{A}, Y_i)}{\text{LNW}(\Omega^{B}, X_i)}} \\

 Like prior work, I assume that operations in semantic vector spaces are reasonably analogous to the targeted cultural-associative processes and that standard embedding models are capable (in principle) of estimating vector spaces that adequately support descriptions of these associations.
 
 # Bias: link between covariance ratios.
# There is overlap where the gender (pink/red) and race (yellow/teal) categories are repeated
fblm_results %>% ggplot(aes(x=as.numeric(b2.cr), y=as.numeric(b3.cr), color=comparison, label=institution)) + geom_text()

An alternative distance-oriented approach is to employ set-valued distance metrics such as the Hausdorff distance, which measures the furthest nearest neighbor in set B of every point in set A (or vice-versa). This is an asymmetric quantity ("the" Hausdorff distance is the maximum of both of these distances) so examining both sides is usually preferable.