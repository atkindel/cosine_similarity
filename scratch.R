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


