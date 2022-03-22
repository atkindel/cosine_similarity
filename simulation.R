library(tidyverse)
library(textdata)
library(gridExtra)

theme_set(theme_bw())


## Construct embeddings

glove_d <- embedding_glove6b(dimensions=200)
glove <- glove_d %>% pivot_longer(contains("d"), names_to="dim")
glove_tokens <- unique(glove$token)
glove_d %<>% select(-token) %>% as.matrix()
rownames(glove_d) <- glove_tokens

# WEAT black/white names and good/bad terms
weat_white <- c("adam", "harry", "josh", "roger", "alan", "frank", "justin", "ryan", "andrew", "jack",
                "matthew", "stephen", "brad", "greg", "paul", "jonathan", "peter", "amanda", "courtney",
                "heather", "melanie", "katie", "betsy", "kristin", "nancy", "stephanie", "ellen", "lauren",
                "colleen", "emily", "megan", "rachel")
weat_black <- c("alonzo", "jamel", "theo", "alphonse", "jerome", "leroy", "torrance", "darnell", "lamar",
                "lionel", "tyree", "deion", "lamont", "malik", "terrence", "tyrone", "lavon", "marcellus",
                "wardell", "nichelle", "shereen", "ebony", "latisha", "shaniqua", "jasmine", "tanisha",
                "tia", "lakisha", "latoya", "yolanda", "malika", "yvette")
weat_good <- c("caress", "freedom", "health", "love", "peace", "cheer", "friend", "heaven", "loyal", "pleasure",
               "diamond", "gentle", "honest", "lucky", "rainbow", "diploma", "gift", "honor", "miracle", "sunrise",
               "family", "happy", "laughter", "paradise", "vacation")
weat_bad <- c("abuse", "crash", "filth", "murder", "sickness", "accident", "death", "grief", "poison", "stink",
              "assault", "disaster", "hatred", "pollute", "tragedy", "bomb", "divorce", "jail", "poverty", "ugly",
              "cancer", "evil", "kill", "rotten", "vomit")

# Get GloVe matrices
weat_white_glove <- scale(glove_d[c(weat_white),])
weat_black_glove <- scale(glove_d[c(weat_black),])
weat_good_glove <- scale(glove_d[c(weat_good),])
weat_bad_glove <- scale(glove_d[c(weat_bad),])

# Look at pairwise cosine distributions for each space over a few subsamples of random size

# Function to compute cosine similarity over a subset of WEAT terms
repl_subs <- function(word_list, ns) {
  glss <- glove_d[word_list,]
  glss <- glss[sample(1:ns),]
  glss <- scale(glss)
  data.frame(lsa::cosine(t(glss))) %>%
    mutate(t1 = rownames(.)) %>%
    pivot_longer(1:last_col(1), names_to="t2") %>%
    filter(t1 != t2) ->  # self cosine is always 1
    ct
  return(as.matrix(ct))
}

# Get one run
n_samp <- 25
wbad <- repl_subs(weat_bad, n_samp)
wgood <- repl_subs(weat_good, n_samp)
wwhite <- repl_subs(weat_white, n_samp)
wblack <- repl_subs(weat_black, n_samp)

rbind(data.frame(wbad) %>% mutate(group="bad"), data.frame(wgood) %>% mutate(group="good"),
      data.frame(wwhite) %>% mutate(group="white"), data.frame(wblack) %>% mutate(group="black")) -> wbg

wbg %>% ggplot(aes(x=as.numeric(value), fill=group, color=group)) + geom_density(alpha=0.25)



# Inner product of corresponding GloVe vectors
gdot <- function(tx, ty) {
  geometry::dot(glove_d[tx,], glove_d[ty,])
}

# Get the cosine similarity distribution for each group comparison
weat_white_glove <- glove_d[c(weat_white),]
weat_black_glove <- glove_d[c(weat_black),]
weat_good_glove <- glove_d[c(weat_good),]
weat_bad_glove <- glove_d[c(weat_bad),]

rbind(data.frame(weat_bad_glove) %>% slice_sample(n=n_samp) %>% mutate(group="bad"),
      data.frame(weat_good_glove) %>% slice_sample(n=n_samp) %>% mutate(group="good"),
      data.frame(weat_white_glove) %>% slice_sample(n=n_samp) %>% mutate(group="white"),
      data.frame(weat_black_glove) %>% slice_sample(n=n_samp) %>% mutate(group="black")) -> wbgb

wbgb_c <- lsa::cosine(t(as.matrix(wbgb %>% select(-group))))
wbgb_d <- cbind(data.frame(wbgb_c), wbgb %>% select(tgroup=group)) %>%
  mutate(term = rownames(.)) %>%
  pivot_longer(1:last_col(2)) %>%
  filter(term != name) %>%
  left_join(wbgb %>% select(ngroup=group) %>% mutate(name=rownames(.)), by="name") %>%
  rowwise() %>%
  mutate(ip = gdot(term, name))


# wbgb_d %>%
#   group_by(tgroup, ngroup) %>%
#   mutate(mv = mean(value)) %>%
#   ggplot(aes(x=value, color=ngroup, fill=ngroup)) +
#   geom_density(alpha=0.25) +
#   geom_rug(alpha=0.5) +
#   geom_vline(aes(xintercept=mv, group=ngroup, color=ngroup), linetype="dashed") +
#   facet_wrap(~tgroup)

# Above plot but with marginal boxplots
p_cos <- function(tdf, tg, inner=F) {
  tdf %>%
    mutate(hv = ifelse(inner, ip, value)) %>%
    group_by(tgroup, ngroup) %>%
    mutate(mv = mean(hv)) %>%
    ungroup() %>%
    filter(tgroup==tg) %>%
    ggplot(aes(x=hv, color=ngroup, fill=ngroup)) +
    geom_density(alpha=0.25) +
    geom_point(aes(y=0.1), alpha=0) +
    geom_vline(aes(xintercept=mv, group=ngroup, color=ngroup), linetype="dashed") +
    labs(x=ifelse(inner, "Inner product", TeX("$cos(\\theta)$")),
         title=paste0("Term group: ", tg)) -> p1
  p1m <- ggMarginal(p1, margin="x", type="boxplot", size=5, groupColour=T, groupFill=T)
  return(p1m)
}

p_cos_all <- function(tdf, inner=F) {
  grid.arrange(p_cos(tdf, "good", inner), p_cos(tdf, "bad", inner),
               p_cos(tdf, "white", inner), p_cos(tdf, "black", inner), ncol=2)
}

p_cos_all(wbgb_d, inner=F)
p_cos_all(wbgb_d, inner=T)

# What happens if you scale the subspaces first?
weat_white_glove_sc <- scale2(glove_d[c(weat_white),])
weat_black_glove_sc <- scale2(glove_d[c(weat_black),])
weat_good_glove_sc <- scale2(glove_d[c(weat_good),])
weat_bad_glove_sc <- scale2(glove_d[c(weat_bad),])

rbind(data.frame(weat_bad_glove_sc) %>% slice_sample(n=n_samp) %>% mutate(group="bad"),
      data.frame(weat_good_glove_sc) %>% slice_sample(n=n_samp) %>% mutate(group="good"),
      data.frame(weat_white_glove_sc) %>% slice_sample(n=n_samp) %>% mutate(group="white"),
      data.frame(weat_black_glove_sc) %>% slice_sample(n=n_samp) %>% mutate(group="black")) -> wbgb_s

wbgb_sc <- lsa::cosine(t(as.matrix(wbgb_s %>% select(-group))))
wbgb_sd <- cbind(data.frame(wbgb_sc), wbgb %>% select(tgroup=group)) %>%
  mutate(term = rownames(.)) %>%
  pivot_longer(1:last_col(2)) %>%
  filter(term != name) %>%
  left_join(wbgb_s %>% select(ngroup=group) %>% mutate(name=rownames(.)), by="name") %>%
  rowwise() %>%
  mutate(ip = gdot(term, name))

p_cos_all(wbgb_sd, inner=F)
p_cos_all(wbgb_sd, inner=T)


# Compare joint distributions
# In the globally scaled distribution it's basically the same line (but heteroskedastic, I am not letting this go lol)
# In the locally scaled distribution it's still heteroskedastic but also there is super clear frequency bias I think?
wbgb_sd %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_jitter(alpha=0.2) +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product")

wbgb_d %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_jitter(alpha=0.2) +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product")

# Use glove index as a frequency rank statistic
frank <- function(glist) {
  data.frame(gidx = which(glove_tokens %in% glist),
             token = glove_tokens[which(glove_tokens %in% glist)]) %>%
    mutate(fr=rank(gidx))
}
rbind(frank(weat_bad),
      frank(weat_good),
      frank(weat_black),
      frank(weat_white)) ->
  frankd

# Alternatively, use the magnitude (Euclidean norm) of the vector as the statistic
megan <- function(glist) {
  data.frame(mag = sapply(1:nrow(glove_d[glist,]), function(x) norm(glove_d[glist,x], type="2")),
             token = glove_tokens[which(glove_tokens %in% glist)])
}

rbind(megan(weat_bad),
      megan(weat_good),
      megan(weat_black),
      megan(weat_white)) ->
  megand

wbgb_sd %>%
  left_join(frankd %>% rename(t_grank = gidx, tfr=fr), by=c("term" = "token")) %>%
  left_join(frankd %>% rename(n_grank = gidx, nfr=fr), by=c("name" = "token")) %>%
  mutate(csize = t_grank / n_grank,
         lcsize = log(csize)) ->
  wbgb_sdr

wbgb_d %>%
  left_join(frankd %>% rename(t_grank = gidx, tfr=fr), by=c("term" = "token")) %>%
  left_join(frankd %>% rename(n_grank = gidx, nfr=fr), by=c("name" = "token")) %>%
  left_join(megand %>% rename(tmag = mag), by=c("term"="token")) %>%
  left_join(megand %>% rename(nmag = mag), by=c("name"="token")) %>%
  mutate(csize = t_grank / n_grank,
         lcsize = log(csize),
         msize = tmag / nmag) ->
  wbgb_dr

# Frequency bias
# In the scaled comparison, there is a robust multiplicative interaction between
#  cosine similarity and frequency imbalance.
# summary(lm(ip ~ value, data=wbgb_sdr))
# summary(lm(ip ~ value, data=wbgb_dr))
# 
# summary(lm(csize ~ value, data=wbgb_sdr))
# summary(lm(csize ~ value, data=wbgb_dr))
# 
# summary(lm(ip ~ value * csize, data=wbgb_sdr))
# summary(lm(ip ~ value * csize, data=wbgb_dr))
# 
# m1 <- MASS::rlm(ip ~ value * csize, data=wbgb_sdr, method="MM")
# m2 <- mgcv::gam(ip ~ s(value) + s(csize) + s(value, csize), data=wbgb_sdr)
# 
# m1b <- MASS::rlm(ip ~ value * csize, data=wbgb_dr, method="MM")
# m2b <- mgcv::gam(ip ~ s(value) + s(csize) + s(value, csize), data=wbgb_dr)
# In general the plots below are a better way of showing this

# Show frequency imbalance as point size
# Look at the way the frequency gradient cuts across the distribution
# This pattern does not appear if you don't scale the vectors locally
wbgb_sdr %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_point(aes(size=lcsize), alpha=0.2) +
  geom_smooth(method="lm") +
  scale_size_binned_area(n.breaks=10) +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       size="Frequency rank ratio",
       color="Term group")

# Plot color as frequency imbalance, use shape for group
wbgb_sdr %>%
  ggplot(aes(x=value, y=ip, shape=tgroup)) +
  geom_point(aes(color=abs(lcsize)), alpha=0.2, size=3) +
  geom_smooth(method="lm") +
  scale_color_viridis_c() +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       color="Frequency rank ratio",
       shape="Term group")

# Split by group
wbgb_sdr %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_point(aes(size=lcsize), alpha=0.2) +
  geom_smooth(method="lm") +
  scale_size_binned_area(n.breaks=10) +
  facet_wrap(~tgroup) +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       size="Frequency rank ratio",
       color="Term group")

# Now look at the globally scaled distribution
wbgb_dr %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_point(aes(size=lcsize), alpha=0.2) +
  geom_smooth(method="lm") +
  scale_size_binned_area(n.breaks=10) +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       size="Frequency rank ratio",
       color="Term group")

# Color-shape plot
wbgb_dr %>%
  ggplot(aes(x=value, y=ip, shape=tgroup)) +
  geom_point(aes(color=abs(lcsize)), alpha=0.2, size=3) +
  geom_smooth(method="lm") +
  scale_color_viridis_c() +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       color="Frequency rank ratio",
       shape="Term group")

# Broken out by group
wbgb_dr %>%
  ggplot(aes(x=value, y=ip, color=tgroup)) +
  geom_point(aes(size=lcsize), alpha=0.2) +
  geom_smooth(method="lm") +
  scale_size_binned_area(n.breaks=10) +
  facet_wrap(~tgroup) +
  labs(x=TeX("$cos(\\theta)$"),
       y="Inner product",
       size="Frequency rank ratio",
       color="Term group")


# Now look at the distribution of cosine differences and frequency rank imbalance ratio
# This implicates 4 terms -- it starts to get really mindbreaking here...
# We drop the degenerate comparisons (self-analogies, mirrored self-analogies)
analogies <- as.data.frame(expand(wbgb_dr, pa=nesting(term, name), pb=nesting(term, name)))
analogies <- cbind(analogies[[1]], analogies[[2]])
colnames(analogies) <- c("pa.term", "pa.name", "pb.term", "pb.name")
analogies %>%
  filter(!((pa.term == pb.term) & (pa.name == pb.name))) %>%
  filter(!((pa.term == pb.name) & (pa.name == pb.term))) %>%
  left_join(wbgb_dr %>%
              select(pa.term=term, pa.name=name,
                     pa.tfr=tfr, pa.nfr=nfr,
                     pa.tmag=tmag, pa.nmag=nmag,
                     pa.trank=t_grank, pa.nrank=n_grank,
                     pa.tgroup=tgroup, pa.ngroup=ngroup,
                     pa.value=value, pa.csize=csize, pa.msize=msize), by=c("pa.term", "pa.name")) %>%
  left_join(wbgb_dr %>%
              select(pb.term=term, pb.name=name,
                     pb.tfr=tfr, pb.nfr=nfr,
                     pb.tmag=tmag, pb.nmag=nmag,
                     pb.trank=t_grank, pb.nrank=n_grank,
                     pb.tgroup=tgroup, pb.ngroup=ngroup,
                     pb.value=value, pb.csize=csize, pb.msize=msize), by=c("pb.term", "pb.name")) %>%
  mutate(vdiff = pa.value - pb.value,
         csize_r = pa.csize/pb.csize,
         msize_r = pa.msize/pb.msize,
         n_rr = pa.nrank / pb.nrank,
         n_mm = pa.nmag / pb.nmag) ->
  analogy_space

headspace <- head(analogy_space, 1000)

# Subset to only the comparisons we care about
# This is s(w, A, B) -- w can be good or bad words, but they are always the same
# It also only compares ordered pairs of names and ordered pairs of good/bad association words
# Note the statistic does not take the inner geometry of the good and bad words into account
# The comparisons we have selected are a tiny fraction of the total space of comparisons
# s(w, A, B) is correlated with the relative frequency ratio in this subspace
analogy_space %>%
  filter(pa.tgroup %in% c("good", "bad") & pa.term == pb.term) %>%  # Good-good and bad-bad comparisons
  filter(pa.ngroup == "white" & pb.ngroup == "black" | pa.ngroup == "black" & pb.ngroup == "white") %>%  # White-black comparisons
  filter(pa.nfr == pb.nfr & pa.tfr == pb.tfr) %>%  # Name-name and attribute-attribute within-rank comparisons
  filter(pa.ngroup == "white") ->  # Absolute difference is one-sided, there is no s(w, B, A) term
  weat_space

# Plot difference function with relative frequency rank ratio as color
# The evenness of the distribution of the comparison throughout the function varies
#  systematically with relative frequency!
weat_space %>%
  # ggplot(aes(x=pa.value, y=pb.value, color=log(csize_r))) +
  ggplot(aes(x=pa.value, y=pb.value, color=msize_r)) +
  geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c()

# Split by term
weat_space %>%
  ggplot(aes(x=pa.value, y=pb.value, color=msize_r)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  facet_wrap(~pa.term) +
  scale_color_viridis_c()

# More uneven comparisons have bigger residuals in the relationship between
#  cosine difference and relative frequency ratio. We can observe a positive
#  relationship between the relative frequency ratio and the cosine difference.
weat_space %>%
  mutate(size_r = resid(lm(vdiff ~ msize_r, data=weat_space))) %>%
  ggplot(aes(x=msize_r, y=vdiff, color=abs(size_r))) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  scale_color_viridis_c()


# Now summing up to S(X, Y, A, B)
weat_space %>%
  group_by(pa.term) %>%
  summarize(swab = mean(vdiff),
            pa.tfr=first(pa.tfr),
            pa.tgroup=first(pa.tgroup)) %>%
  group_by(pa.tgroup) %>%
  mutate(psum = sum(swab)) ->
  weat_summ

# Note that the sum upweights the larger differences
weat_summ %>%
  left_join(weat_summ, by=c("pa.tfr")) %>%
  filter(pa.tgroup.x != pa.tgroup.y) %>%
  filter(psum.y > 3) %>%
  left_join(frankd %>% mutate(pa.term.x.fd = gidx), by=c("pa.term.x" = "token")) %>%
  left_join(frankd %>% mutate(pa.term.y.fd = gidx), by=c("pa.term.y" = "token")) %>%
  mutate(t_rr = pa.term.x.fd / pa.term.y.fd) ->
  weat_summ_2

# Plot terms of difference function on response surface
# Color term-term differences in association by their relative frequency
weat_summ_2 %>%
  ggplot(aes(y=swab.x, x=swab.y, color=t_rr, label=paste0(pa.term.x, " : ", pa.term.y))) +
  geom_text(size=4) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c() +
  labs(x="good", y="bad")

# Plot the relationship between relative frequency and the WEAT statistic components
weat_summ_2 %>%
  ggplot(aes(y=swab.x - swab.y, x=t_rr, label=paste0(pa.term.x, " : ", pa.term.y))) +
  geom_text(size=4) +
  geom_smooth() +
  scale_color_viridis_c()



# Next, look at the geometry of the comparison vectors
# The cosine similarity is the column 

# Get the cosine similarity distribution for each group comparison
weat_white_glove <- glove_d[c(weat_white),]
weat_black_glove <- glove_d[c(weat_black),]
weat_good_glove <- glove_d[c(weat_good),]
weat_bad_glove <- glove_d[c(weat_bad),]

# OPTIONALLY, reorder the rows
weat_white_b <- sample(c(weat_white))
weat_black_b <- sample(c(weat_black))
weat_good_b <- sample(c(weat_good))
weat_bad_b <- sample(c(weat_bad))
weat_white_glove_b <- glove_d[c(weat_white_b),]
weat_black_glove_b <- glove_d[c(weat_black_b),]
weat_good_glove_b <- glove_d[c(weat_good_b),]
weat_bad_glove_b <- glove_d[c(weat_bad_b),]

get_similarity_loading <- function(mat_a, mat_b, ndim=25, zt=0.0000001, operator="-") {
  tdf <- NULL
  if(operator == "+") {
    tdf <- data.frame(mat_a[1:ndim,] + mat_b[1:ndim,])
  } else {
    tdf <- data.frame(mat_a[1:ndim,] - mat_b[1:ndim,])
  }
  tdf %>%
    mutate(nom = rownames(.)) %>%
    pivot_longer(1:last_col(1)) %>%
    mutate(value2=ifelse(abs(value) <= zt, NA, value)) ->  # Mark effective zeros as missing
    mminus
  levels(mminus$name) <- unique(mminus$name)
  return(mminus)
}

plot_similarity_loading <- function(dm, zero_missing=F) {
  dm %>%
    rowwise() %>%
    mutate(value = ifelse(zero_missing, value2, value)) %>%
    ggplot(aes(x=nom, y=name, fill=abs(value))) +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_y_discrete(labels=unique(wbg_minus$name))
}

wbg_minus <- get_similarity_loading(weat_black_glove, weat_good_glove)
wbg_plus <- get_similarity_loading(weat_black_glove, weat_good_glove, operator="+")

# Pretty clear that these comparisons target different dimensions in general
plot_similarity_loading(wbg_minus)
plot_similarity_loading(wbg_plus)

# wbg_minus %>%
#   group_by(name) %>%
#   mutate(dim_norm = norm(value, "2"),
#          val_ds = value/dim_norm) %>%
#   group_by(nom) %>%
#   mutate(term_norm = norm(value, "2"),
#          val_ts = value/term_norm) ->
#   wbg_minus_marg

# weat_white_glove[1:25,] - weat_black_glove[1:25,] -> test_a -> test
test <- test_b <- t(sapply(1:25, function(x) scale(weat_white_glove[x,]) - scale(weat_black_glove[x,])))
rownames(test) <- paste0(rownames(weat_white_glove[1:25,]), "_", rownames(weat_black_glove[1:25,]))
test_good <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], test[x,]))
test_bad <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], test[x,]))
dnorm <- apply(test, 1, norm, "2")  # spectral norm of difference vector (largest singular value)
good_black <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_black_glove[x,]))
bad_black <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], weat_black_glove[x,]))
good_white <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_white_glove[x,]))
bad_white <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], weat_white_glove[x,]))
black_white <- sapply(1:25, function(x) lsa::cosine(weat_black_glove[x,], weat_white_glove[x,]))
good_bad <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_bad_glove[x,]))

# Function to get angle between D and W in linear subspace spanned by AB
proj_angle <- function(i, j=i, wa, wb, ww) {
  # Find subspace in Rp
  data.frame(a=wa[i,], b=wb[i,]) -> atest
  svd(atest)$u -> wbb
  
  # Project W onto AB
  wbb_proj <- wbb %*% solve(t(wbb) %*% wbb) %*% t(wbb)
  wbb_proj_ww <- as.vector(wbb_proj %*% ww[j,])

  # Compute angle
  # lsa::cosine(wbb_proj_ww, ww[i,])
  lsa::cosine(as.vector(scale(wa[i,]) - scale(wb[i,])), wbb_proj_ww)
}

wba_good <- sapply(1:25, proj_angle, wa=weat_white_glove, wb=weat_black_glove, ww=weat_good_glove)
wba_bad <- sapply(1:25, proj_angle, wa=weat_white_glove, wb=weat_black_glove, ww=weat_bad_glove)

cbind.data.frame(test_good, test_bad, good_black, bad_black, good_white, bad_white, black_white, good_bad, dnorm, wba_good, wba_bad) %>%
  mutate(tgn = test_good * dnorm, good_diff = good_white - good_black,
         tbn = test_bad * dnorm, bad_diff = bad_white - bad_black) ->
  test_results

# Let the sentiment words be X=good and Y=bad, and the attribute words be A=white and B=black
# The difference in cosine similarity should be the cosine similarity of the
#  difference times a known constant ||D||, but it isn't.
# It's noticeably more off the larger the standard basis difference.
# This is because the x-axis cosines are not in the same subspace as the y-axis cosine
test_results %>%
  ggplot(aes(x=good_diff, y=tgn)) +
  geom_point(aes(size=abs(wba_good))) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
       y=TeX("$||D||cos(\\theta_{X_i, A_i - B_i})$"))

test_results %>%
  ggplot(aes(x=bad_diff, y=tbn, size=wba_bad)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta_{Y_i, A_i}) - cos(\\theta_{Y_i, B_i})$"),
       y=TeX("$||D||cos(\\theta_{Y_i, A_i - B_i})$"))

# Get angle between linear subspace AB and word vector W
# data.frame(white=weat_white_glove[1,], black=weat_black_glove[1,]) -> test
# svd(test)$u -> wbb
# wbb_proj <- wbb %*% solve(t(wbb) %*% wbb) %*% t(wbb)
# wbb_proj_good <- as.vector(wbb_proj %*% weat_good_glove[1,])
# wbb_proj_bad <- as.vector(wbb_proj %*% weat_bad_glove[1,])
# lsa::cosine(wbb_proj_good, weat_good_glove[1,])



# weat_white_glove[1:25,] - weat_black_glove[1:25,] -> test_a -> test
test <- test_b <- t(sapply(1:25, function(x) scale(weat_white_glove[x,]) - scale(weat_black_glove[x,])))
rownames(test) <- paste0(rownames(weat_white_glove[1:25,]), "_", rownames(weat_black_glove[1:25,]))
test_good <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], test[x,]))
test_bad <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], test[x,]))
dnorm <- apply(test, 1, norm, "2")  # spectral norm of difference vector (largest singular value)
good_black <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_black_glove[x,]))
bad_black <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], weat_black_glove[x,]))
good_white <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_white_glove[x,]))
bad_white <- sapply(1:25, function(x) lsa::cosine(weat_bad_glove[x,], weat_white_glove[x,]))
black_white <- sapply(1:25, function(x) lsa::cosine(weat_black_glove[x,], weat_white_glove[x,]))
good_bad <- sapply(1:25, function(x) lsa::cosine(weat_good_glove[x,], weat_bad_glove[x,]))


# Bilinear form
# Difference of dot products:
sum(weat_white_glove[1,]/norm(weat_white_glove[1,], "2") *
      weat_good_glove[1,]/norm(weat_good_glove[1,], "2")) -
  sum(weat_black_glove[1,]/norm(weat_black_glove[1,], "2") *
        weat_good_glove[1,]/norm(weat_good_glove[1,], "2"))
# Equals dot product of differences:
sum(weat_good_glove[1,]/norm(weat_good_glove[1,], "2") *
      (weat_white_glove[1,]/norm(weat_white_glove[1,], "2") - weat_black_glove[1,]/norm(weat_black_glove[1,], "2")))

scale2 <- function(vc) {
  return(vc / norm(vc, "2"))
}

# Now look at all of the angle projections that go into the statistic
# The point of these is to show that cos(WA) + cos(WB) = ||D||cos(WD) where D=A/||A|| + B/||B||
# The scale dependence of the sum can be thought about in terms of the
#  scalar projection of this normalized sum vector onto the comparison word W
expand_grid(i=1:32, j=1:25) %>%
  rowwise() %>%
  mutate(wname = weat_white[i],
         bname = weat_black[i],
         wterm = weat_good[j],
         wn_norm = norm(weat_white_glove[i,], "2"),
         bn_norm = norm(weat_black_glove[i,], "2"),
         wt_norm = norm(weat_good_glove[j,], "2"),
         wbg_adiff = lsa::cosine(weat_white_glove[i,], weat_good_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_good_glove[j,]),
         wbg_pang = proj_angle(i, j, wa=weat_white_glove, wb=weat_black_glove, ww=weat_good_glove),
         wbg_wbang = lsa::cosine(weat_black_glove[i,], weat_white_glove[i,]),
         wbg_wbsang = lsa::cosine(scale2(weat_black_glove[i,]), scale2(weat_white_glove[i,])),
         wbg_csang = lsa::cosine(as.vector(scale2(weat_white_glove[i,]) + scale2(weat_black_glove[i,])), weat_good_glove[j,]),
         wbg_cuang = lsa::cosine(weat_white_glove[i,] - weat_black_glove[i,], weat_good_glove[j,]),
         wbg_csdnorm = norm(as.vector(scale2(weat_white_glove[i,]) + scale2(weat_black_glove[i,])), "2"),
         wbg_cmdnorm = norm(scale2(weat_white_glove[i,] - weat_black_glove[i,]), "2"),
         csmd_ratio = wbg_csdnorm / wbg_cmdnorm,
         wbg_cudnorm = norm(weat_white_glove[i,] - weat_black_glove[i,], "2"),
         wbg_dnang = 1 - (wbg_csdnorm ^ 2)/2,
         wbg_angbynorm = wbg_csang * wbg_csdnorm,
         pairwise_3_sum = lsa::cosine(weat_white_glove[i,], weat_good_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_good_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_white_glove[i,]),
         pairwise_3s_bilin = wbg_csdnorm * wbg_csang + wbg_wbang
         # sv_norm_ratio = svd(rbind(weat_black_glove[i,],weat_white_glove[i,]))$d[1] / norm(weat_black_glove[i,] + weat_white_glove[i,], "2"),
         # singular_ratio = svd(rbind(weat_black_glove[i,],weat_white_glove[i,]))$d[1] /svd(rbind(weat_black_glove[i,], weat_white_glove[i,]))$d[2],
         # norm_ratio = norm(scale2(weat_black_glove[i,]) + scale2(weat_white_glove[i,]), "2") / norm(scale2(weat_black_glove[i,]) - scale2(weat_white_glove[i,]), "2"),
         # sv_l1_align = lsa::cosine(svd(rbind(weat_black_glove[i,], weat_white_glove[i,]))$v[,1], weat_black_glove[i,] + weat_white_glove[i,]),
         # vmax_l1_align = lsa::cosine((svd(rbind(weat_black_glove[1,], weat_white_glove[1,]))$v %*% varimax(svd(rbind(weat_black_glove[1,], weat_white_glove[1,]))$v)$rotmat)[,1], weat_black_glove[i,] + weat_white_glove[i,])
         ) ->
  angle_projections_good

expand_grid(i=1:32, j=1:25) %>%
  rowwise() %>%
  mutate(wname = weat_white[i],
         bname = weat_black[i],
         wterm = weat_bad[j],
         wn_norm = norm(weat_white_glove[i,], "2"),
         bn_norm = norm(weat_black_glove[i,], "2"),
         wt_norm = norm(weat_bad_glove[j,], "2"),
         wbg_adiff = lsa::cosine(weat_white_glove[i,], weat_bad_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_bad_glove[j,]),  # cos(WA) + cos(WB)
         wbg_pang = proj_angle(i, j, wa=weat_white_glove, wb=weat_black_glove, ww=weat_bad_glove),
         wbg_wbang = lsa::cosine(weat_black_glove[i,], weat_white_glove[i,]),
         wbg_wbsang = lsa::cosine(scale2(weat_black_glove[i,]), scale2(weat_white_glove[i,])),
         wbg_csang = lsa::cosine(as.vector(scale2(weat_white_glove[i,]) + scale2(weat_black_glove[i,])), weat_bad_glove[j,]),  # cos(WD)
         wbg_cuang = lsa::cosine(weat_white_glove[i,] - weat_black_glove[i,], weat_bad_glove[j,]),
         wbg_csdnorm = norm(as.vector(scale2(weat_white_glove[i,]) + scale2(weat_black_glove[i,])), "2"),  # ||D||
         wbg_cmdnorm = norm(scale2(weat_white_glove[i,] - weat_black_glove[i,]), "2"),
         csmd_ratio = wbg_csdnorm / wbg_cmdnorm,
         wbg_cudnorm = norm(weat_white_glove[i,] - weat_black_glove[i,], "2"),
         wbg_dnang = 1 - (wbg_csdnorm ^ 2)/2,
         wbg_angbynorm = wbg_csang * wbg_csdnorm,
         pairwise_3_sum = lsa::cosine(weat_white_glove[i,], weat_bad_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_bad_glove[j,]) + lsa::cosine(weat_black_glove[i,], weat_white_glove[i,])) ->
  angle_projections_bad

# Compare difference to angle -- clearly off (even if you scale the D vectors)
# This is easy to screw up!
angle_projections_good %>%
  mutate(namepair = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=wbg_adiff, y=wbg_csang, label=wterm)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
       y=TeX("$cos(\\theta_{X_i, A_i - B_i})$"))

# Scaled separately (this is how the comparison is constructed)
# This should be right on (it is)
angle_projections_good %>%
  mutate(namepair = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_csdnorm, label=wterm)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x=latex2exp::TeX("$cos(\\theta_{X_i, A_i}) + cos(\\theta_{X_i, B_i})$"),
       y=latex2exp::TeX("$||D||cos(\\theta_{X_i, A_i + B_i})$"))

# Same story for the bad terms
angle_projections_bad %>%
  mutate(namepair = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_csdnorm, label=wterm)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x=latex2exp::TeX("$cos(\\theta_{Y_i, A_i}) - cos(\\theta_{Y_i, B_i})$"),
       y=latex2exp::TeX("$||D||cos(\\theta_{Y_i, A_i - B_i})$"))

# # Scaled jointly
# angle_projections_good %>%
#   mutate(namepair = paste0(wname, " - ", bname)) %>%
#   ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_cmdnorm, label=wterm)) +
#   geom_point() +
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   geom_smooth(method="lm") +
#   labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
#        y=TeX("$||D_{\\psi}||cos(\\theta_{X_i, A_i - B_i})$"))

# # Unscaled
# angle_projections_good %>%
#   mutate(namepair = paste0(wname, " - ", bname)) %>%
#   ggplot(aes(x=wbg_adiff, y=wbg_cuang * wbg_cudnorm, label=wterm)) +
#   geom_point() +
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   geom_smooth(method="lm") +
#   labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
#        y=TeX("$cos(\\theta_{X_i, A_i - B_i})$"),
#        title="Unscaled")

# # Both of the above plots show that neither Dnorm alone works; the key is their ratio
# # NOTE: This is the killer figure!
# angle_projections_good %>%
#   mutate(namepair = paste0(wname, " - ", bname)) %>%
#   ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_csdnorm / wbg_cmdnorm, label=wterm)) +
#   geom_point() +
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   geom_smooth(method="lm") +
#   labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
#        y=TeX("$\\frac{||D_{\\phi}||}{||D_{\\psi}||}cos(\\theta_{X_i, A_i - B_i})$"))

# This holds in every name comparison subspace
angle_projections_good %>%
  mutate(namepair = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_csdnorm, label=wterm)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_point() +
  labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
       y=TeX("$||D||cos(\\theta_{X_i, A_i - B_i})$")) +
  facet_wrap(~namepair)

# Bad terms by pair subspaces
angle_projections_bad %>%
  mutate(namepair = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=wbg_adiff, y=wbg_csang * wbg_csdnorm, label=wterm)) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_point() +
  labs(x=TeX("$cos(\\theta_{Y_i, A_i}) - cos(\\theta_{Y_i, B_i})$"),
       y=TeX("$||D||cos(\\theta_{Y_i, A_i - B_i})$")) +
  facet_wrap(~namepair)

# Look at distribution of ||D||
# (Does not depend on choice of W.)
angle_projections_good %>%
  ggplot(aes(x=wbg_csdnorm)) +
  geom_density()

# Look at ||D|| as a function of ||A|| and ||B||
angle_projections_good %>%
  ggplot(aes(x=wn_norm, y=bn_norm, color=wbg_csdnorm)) +
  geom_point(size=3) +
  geom_smooth(method="lm")

# Look at the cosine similarity cosWD as a function of ||A|| and ||B||
# This has two levels: joint distribution over ||A||, ||B||; distribution within each "point"
angle_projections_good %>%
  ggplot(aes(x=wn_norm, y=bn_norm, color=wbg_csang)) +
  geom_point(size=3) +
  scale_color_viridis_c()

# The underlying distributions are generally different
angle_projections_good %>%
  mutate(pairn = paste0(wname, ", ", bname)) %>%
  ggplot(aes(x=wbg_csang, group=pairn, fill=pairn)) +
  geom_density(alpha=0.1)

# The three-way paired sum is also ultra weird
# cosAB + cosBC + cosAC
# There is a cloud of these norm triples; the distribution of 
angle_projections_good %>%
  ggplot(aes(x=wn_norm, y=bn_norm)) +
  geom_point(aes(color=pairwise_3_sum, size=pairwise_3_sum, alpha=wt_norm)) +
  facet_wrap(~wterm) +
  labs(color="Pairwise sum (k=3)", size="") +
  theme(legend.position="bottom")

# Each comparison has a different distribution
angle_projections_good %>%
  ggplot(aes(x=pairwise_3_sum, color=wt_norm)) +
  geom_density() +
  facet_wrap(~wterm) +
  labs(color="||W||") +
  theme(legend.position="bottom")




# That was the distribution of s(w, A, B)
# Now look at the full four-way comparison set
# Think about this in terms of the three-angle coordinate system in the paper
i_names <- length(weat_black)
j_terms <- length(weat_good)
expand_grid(i=1:i_names, j=1:j_terms) %>%
  rowwise() %>%
  mutate(wname = weat_white_b[i],
         bname = weat_black_b[i],
         gterm = weat_good_b[j],
         bterm = weat_bad_b[j],
         wn_norm = norm(weat_white_glove[i,], "2"),
         bn_norm = norm(weat_black_glove[i,], "2"),
         wt_norm = norm(weat_good_glove[j,], "2"),
         bt_norm = norm(weat_bad_glove[j,], "2"),
         ab_dnorm = norm(as.vector(scale2(weat_white_glove_b[i,]) - scale2(weat_black_glove_b[i,])), "2"),
         xy_tnorm = norm(as.vector(scale2(weat_good_glove_b[j,]) - scale2(weat_bad_glove_b[j,])), "2"),
         td_csim = lsa::cosine(as.vector(scale2(weat_white_glove_b[i,]) - scale2(weat_black_glove_b[i,])),
                               as.vector(scale2(weat_good_glove_b[j,]) - scale2(weat_bad_glove_b[j,]))),
         sabxy_fact = ab_dnorm * xy_tnorm * td_csim,
         ab_csim = 1 - (ab_dnorm ^ 2)/2,
         xy_csim = 1 - (xy_tnorm ^ 2)/2,
         scaling_ratio = ab_dnorm * xy_tnorm,
         abxy_csratio = ab_csim / xy_csim,
         paired_sum=lsa::cosine(weat_white_glove[i,], weat_black_glove[i,]) + lsa::cosine(weat_good_glove[j,], weat_bad_glove[j,]),
         wb_cs = lsa::cosine(weat_white_glove[i,], weat_black_glove[i,]),
         gb_cs = lsa::cosine(weat_good_glove[j,], weat_bad_glove[j,])) ->
  analogy_projections
analogy_projections %<>% group_by(i) %>% mutate(im=mean(sabxy_fact)) %>% group_by(j) %>% mutate(jm = mean(sabxy_fact))

# Look at paired sum as a function of the ratio of the norm products
# Also look at the corresponding cosine similarities
analogy_projections %>%
  ggplot(aes(x=(bt_norm*wt_norm) / (wn_norm*bn_norm),
             y=paired_sum,
             color=wb_cs, size=gb_cs)) +
  geom_point()

# The comparison does not have equal variance across pairs of pairs (next 2 plots)
analogy_projections %>%
  mutate(nd = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=sabxy_fact, fill=nd)) +
  geom_density(alpha=0.25) +
  geom_histogram(alpha=0.8) +
  geom_vline(aes(xintercept=mean(sabxy_fact))) + 
  geom_vline(aes(xintercept=im, color=nd), linetype="dashed") +
  facet_wrap(~nd) +
  theme(legend.position="none") +
  labs(x=TeX("$\\phi(cos(\\theta_{AB}))\\cdot\\phi(cos(\\theta_{XY}))\\cdot cos(\\theta_{TD})$"))

analogy_projections %>%
  mutate(nt = paste0(gterm, " - ", bterm)) %>%
  ggplot(aes(x=sabxy_fact, fill=nt)) +
  geom_density(alpha=0.25) +
  geom_histogram(alpha=0.8) +
  geom_vline(aes(xintercept=mean(sabxy_fact))) + 
  geom_vline(aes(xintercept=jm, color=nt), linetype="dashed") +
  facet_wrap(~nt) +
  theme(legend.position="none") +
  labs(x=TeX("$\\phi(cos(\\theta_{AB}))\\cdot\\phi(cos(\\theta_{XY}))\\cdot cos(\\theta_{TD})$"))

# The statistic is a sum of means of this product:
analogy_projections %>%
  ggplot(aes(x=ab_dnorm, y=xy_tnorm, color=td_csim)) +
  geom_point(size=3) +
  scale_color_viridis_c()

# Smoothed
analogy_projections %>%
  ggplot(aes(x=ab_dnorm, y=xy_tnorm)) +
  geom_density_2d_filled(bins=80) +
  theme(legend.position = "none")

# Another way of seeing it
analogy_projections %>%
  mutate(nd = paste0(wname, " - ", bname)) %>%
  ggplot(aes(x=sabxy_fact, fill=nd)) +
  geom_density(alpha=0.25) +
  geom_vline(aes(xintercept=mean(sabxy_fact))) + theme(legend.position="none")

# And another
# When the difference vectors are more disaligned, the scale is systematically larger!
# These are actually clustered by comparison too (as the norms show)
analogy_projections %>%
  ggplot(aes(x=scaling_ratio, y=td_csim, color=xy_tnorm)) +
  # geom_smooth(method="lm", color="grey80", linetype="dashed") + 
  geom_smooth() +
  geom_point(aes(size=abs(sabxy_fact))) +
  scale_color_viridis_c() +
  labs(x=latex2exp::TeX("$\\phi(cos(\\theta_{AB}))\\cdot\\phi(cos(\\theta_{XY}))$"),
       y=latex2exp::TeX("$cos(\\theta_{TD})$"),
       color=latex2exp::TeX("$||D_i||$"),
       size=latex2exp::TeX("$\\psi_{ij}$"))->
  p1

analogy_projections %>%
  ggplot(aes(x=scaling_ratio, y=td_csim, color=ab_dnorm)) +
  # geom_smooth(method="lm", color="grey80", linetype="dashed") + 
  geom_smooth() +
  geom_point(aes(size=abs(sabxy_fact))) +
  scale_color_viridis_c() +
  labs(x=latex2exp::TeX("$\\phi(cos(\\theta_{AB}))\\cdot\\phi(cos(\\theta_{XY}))$"),
       y=latex2exp::TeX("$cos(\\theta_{TD})$"),
       color=latex2exp::TeX("$||T_j||$"),
       size=latex2exp::TeX("$\\psi_{ij}$"))->
  p2

grid.arrange(p1, p2)

analogy_projections %>%
  ggplot(aes(y=td_csim, x=sabxy_fact, size=scaling_ratio, color=log(abs(abxy_csratio)))) +
  geom_smooth(method="lm") +
  geom_point(alpha=0.5) +
  scale_color_viridis_c()

# The comparison vectors are not normed, so there is information about their relative scale in the statistic
# The alignment we want is td_csim but this is scaled by a ratio of their norms
# We know the exact functional form of the ratio: 4*sin(0.5*theta_ab)*sin(0.5*theta_xy)
# The product of the half-angle sines of the angles theta_ab and theta_xy
scrat_f <- function(cs_ab, cs_xy) {
  return(4 * sqrt(2 - 2*cs_ab) * sqrt(2 - 2*cs_xy))
}

analogy_projections %>%
  mutate(analytic_scale_ratio = scrat_f(ab_csim, xy_csim)) %>%
  ggplot(aes(x=analytic_scale_ratio, y=scaling_ratio)) +
  geom_point()

data.frame(x=seq(-1, 1, by=0.01), y=seq(-1, 1, by=0.01)) %>%
  expand(x, y) %>%
  mutate(scale_rf = scrat_f(x, y)) %>%
  ggplot(aes(x, y, fill=scale_rf)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x=latex2exp::TeX("$\\alpha$"),
       y=latex2exp::TeX("$\\beta$"),
       fill="Weight")


  
  
# Relationship to index/rank/magnitude order effect

# Look at the component products against the ratio of cosine similarities
analogy_projections %>%
  ggplot(aes(x=sabxy_fact, y=abxy_csratio)) +
  geom_point(alpha=0.5)

# Append rank order data
analogy_projections %>%
  left_join(frankd %>% transmute(wname=token, wname_gidx=gidx, wname_rank=fr), by="wname") %>%
  left_join(frankd %>% transmute(bname=token, bname_gidx=gidx, bname_rank=fr), by="bname") %>%
  left_join(frankd %>% transmute(gterm=token, gterm_gidx=gidx, gterm_rank=fr), by="gterm") %>%
  left_join(frankd %>% transmute(bterm=token, bterm_gidx=gidx, bterm_rank=fr), by="bterm") %>%
  mutate(ab_idxratio = wname_gidx/bname_gidx,
         xy_idxratio = gterm_gidx/bterm_gidx,
         ab_rankratio = wname_rank/bname_rank,
         xy_rankratio = gterm_rank/bterm_rank,
         whitenorm = norm(weat_white_glove[i,], "2"),
         blacknorm = norm(weat_black_glove[i,], "2"),
         goodnorm = norm(weat_good_glove[j,], "2"),
         badnorm = norm(weat_bad_glove[j,], "2"),
         ab_mratio = norm(weat_white_glove[i,], "2")/norm(weat_black_glove[i,], "2"),
         xy_mratio = norm(weat_good_glove[j,], "2")/norm(weat_bad_glove[j,], "2"),
         freq_ratio1 = ab_rankratio/xy_rankratio,
         freq_ratio2 = ab_idxratio/xy_idxratio,
         mag_ratio = ab_mratio/xy_mratio) ->
  analogy_rank_projections

# First, just look at the component vector norms against the difference vector norm and the cosine similarity
# Does this look like a well-behaved statistic? Not really
analogy_rank_projections %>%
  ggplot(aes(x=whitenorm, y=blacknorm)) +
  geom_point(aes(color=ab_dnorm, size=ab_csim)) +
  geom_abline(slope=1,intercept=0,linetype="dashed")
analogy_rank_projections %>%
  ggplot(aes(x=whitenorm * blacknorm, y=ab_csim, color=ab_dnorm)) +
  geom_point()

analogy_rank_projections %>%
  ggplot(aes(x=goodnorm, y=badnorm)) +
  geom_point(aes(color=xy_tnorm, size=xy_csim)) +
  geom_abline(slope=1,intercept=0,linetype="dashed")
analogy_rank_projections %>%
  ggplot(aes(x=goodnorm * badnorm, y=xy_csim, color=xy_tnorm)) +
  geom_point()

# The scaling ratio is systematically larger when the comparison is more elliptical
# This plot shows the logratio of index ratios
analogy_rank_projections %>%
  ggplot(aes(x=log(freq_ratio2), y=scaling_ratio)) +
  geom_point() +
  geom_smooth()

# We can also look at the logratio of magnitude ratios of the underlying vectors
analogy_rank_projections %>%
  ggplot(aes(x=log(mag_ratio), y=scaling_ratio)) +
  geom_point() +
  geom_smooth(method="lm")



# Random three-way and four-way comparisons
# Construct independently and identically drawn sets A, B, X, Y from GloVe
# Set a term rarity window and potentially make B and Y more rare on average
# It will be convenient to index A,B and X,Y separately
rarity_imin <- 1
rarity_imax <- 50000
by_offset <- 2000
k <- 100

sA <- sample(glove_tokens[rarity_imin:rarity_imax], k)
sB <- sample(glove_tokens[rarity_imin:rarity_imax], k)
# sX <- sample(glove_tokens[rarity_imin:rarity_imax], k)
# sY <- sample(glove_tokens[rarity_imin:rarity_imax], k)

expand_grid(i=1:k, j=1:k) %>%
  rowwise() %>%
  mutate(aterm = sA[i],
         bterm = sB[i],
         xterm = sX[j],
         yterm = sY[j],
         anorm = norm(glove_d[aterm,], "2"),
         bnorm = norm(glove_d[bterm,], "2"),
         xnorm = norm(glove_d[xterm,], "2"),
         ynorm = norm(glove_d[yterm,], "2"),
         ab_cs = lsa::cosine(glove_d[aterm,], glove_d[bterm,]),
         xy_cs = lsa::cosine(glove_d[xterm,], glove_d[yterm,]),
         ab_dnorm = norm(as.vector(scale2(glove_d[aterm,]) - scale2(glove_d[bterm,])), "2"),
         xy_tnorm = norm(as.vector(scale2(glove_d[xterm,]) - scale2(glove_d[yterm,])), "2"),
         dx_csang = lsa::cosine(as.vector(scale2(glove_d[aterm,]) + scale2(glove_d[bterm,])), glove_d[xterm,]),
         dy_csang = lsa::cosine(as.vector(scale2(glove_d[aterm,]) + scale2(glove_d[bterm,])), glove_d[yterm,]),
         abx_sum = lsa::cosine(glove_d[aterm,], glove_d[xterm,]) + lsa::cosine(glove_d[bterm,], glove_d[xterm,]),
         aby_sum = lsa::cosine(glove_d[aterm,], glove_d[yterm,]) + lsa::cosine(glove_d[bterm,], glove_d[yterm,]),
         inpair_sum = lsa::cosine(glove_d[aterm,], glove_d[xterm,]) + lsa::cosine(glove_d[bterm,], glove_d[xterm,]) + lsa::cosine(glove_d[bterm,], glove_d[aterm,]),
         expair_sum=lsa::cosine(glove_d[aterm,], glove_d[bterm,]) + lsa::cosine(glove_d[xterm,], glove_d[yterm,]),
         twoway_sum=lsa::cosine(glove_d[aterm,], glove_d[xterm,]) + lsa::cosine(glove_d[bterm,], glove_d[xterm,]) + lsa::cosine(glove_d[aterm,], glove_d[yterm,]) + lsa::cosine(glove_d[bterm,], glove_d[yterm,])) ->
  random_projections

random_projections %>%
  ggplot() +
  geom_point(aes(x=anorm, y=bnorm, color=inpair_sum)) +
  scale_color_viridis_c()

random_projections %>%
  ggplot() +
  geom_point(aes(x=xnorm, y=ynorm, color=inpair_sum)) +
  scale_color_viridis_c()



# It can be shown that the sampling distribution of the cosine similarity of a set
# of random vectors is non-normal (it is positively skewed). This has been known since
# the late 19th century and was shown for the Pearson product-moment correlation
# coefficient by Fisher (1928) and subsequently Hotelling (1953).

# Number of vectors in each subspace (assuming equal sized groups)
# When k is lower, the skew is worse; fixed word lists have observably bad geometry
# It's also interesting that the distribution is a function of the amount of data we have.
plot_cs_distributions <- function(k, rmin=1, rmax=nrow(glove_d)) {
  sA <- sample(glove_tokens[rmin:rmax], k)
  sB <- sample(glove_tokens[rmin:rmax], k)
  
  expand_grid(i=1:k, j=1:k) %>%
    rowwise() %>%
    mutate(aterm = sA[i],
           bterm = sB[j],
           biterm = sA[j],
           ab_cs = lsa::cosine(glove_d[aterm,], glove_d[bterm,]),
           abi_cs = lsa::cosine(glove_d[aterm,], glove_d[biterm,])) ->
    random_angles

  # Look at distribution of full graph weight set comparing two vector subspaces
  random_angles %>%
    ggplot(aes(x=ab_cs)) +
    geom_density() +
    geom_vline(aes(xintercept=mean(ab_cs)), linetype="dashed") +
    labs(x="cos(a, b)",
         title="1. Two-way sum") ->
    ra_p1

  # Distribution of cosine similarities taken externally pairwise.
  # A and B are different.
  random_angles %>%
    filter(i == j) %>%
    ggplot(aes(x=ab_cs)) +
    geom_density() +
    geom_vline(aes(xintercept=mean(ab_cs)), linetype="dashed") +
    labs(x="cos(a, b)",
         title="2. Diagonal sum") ->
    ra_p2
  
  # Distribution of cosine similarities taken internally pairwise.
  # A and B are equal.
  random_angles %>%
    filter(i > j) %>%
    ggplot(aes(x=abi_cs)) +
    geom_density() +
    geom_vline(aes(xintercept=mean(abi_cs)), linetype="dashed") +
    labs(x="cos(a, b)",
         title="3. Off-diagonal sum") ->
    ra_p3
  
  grid.arrange(ra_p1, ra_p2, ra_p3, ncol=3)
}

# Plot sum distributions for a range of choices of k
# Larger values take quite a bit longer.
plot_cs_distributions(k=10)
plot_cs_distributions(k=25)
plot_cs_distributions(k=50)
plot_cs_distributions(k=75)
plot_cs_distributions(k=100)
plot_cs_distributions(k=200)







rarity_imin <- 1
rarity_imax <- nrow(glove_d)
k <- 10
# p <- 100  # TODO: What happens when you only use the top d dimensions in the underlying embedding space?

sA <- sample(glove_tokens[rarity_imin:rarity_imax], k)
sB <- sample(glove_tokens[rarity_imin:rarity_imax], k)
expand_grid(i=1:k, j=1:k) %>%
  rowwise() %>%
  mutate(aterm = sA[i],
         bterm = sB[j],
         biterm = sA[j],
         aiterm = sB[i],
         anorm = norm(glove_d[aterm,], "2"),
         bnorm = norm(glove_d[bterm,], "2"),
         ab_cs = lsa::cosine(glove_d[aterm,], glove_d[bterm,]),
         abi_cs = lsa::cosine(glove_d[aterm,], glove_d[biterm,]),
         aib_cs = lsa::cosine(glove_d[aiterm,], glove_d[bterm,]),
         ab_ip = glove_d[aterm,] %*% glove_d[bterm,],
         plus_dnorm = norm(as.vector(scale2(glove_d[aterm,]) + scale2(glove_d[bterm,])), "2"),
         minus_dnorm = norm(as.vector(scale2(glove_d[aterm,]) - scale2(glove_d[bterm,])), "2")) ->
  random_angles

# Polar sine of the original angle of the n-parallelotope implied by each vector subspace
determinant(matrix(random_angles_50$abi_cs, k, k), logarithm=F)
determinant(matrix(random_angles_50$aib_cs, k, k), logarithm=F)

# Look at distribution of full graph weight set comparing two vector subspaces
random_angles %>%
  ggplot(aes(x=ab_cs)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(ab_cs)), linetype="dashed") +
  labs(x="cos(a, b)",
       title="1. Two-way sum") ->
  ra_p1

# # Look from each i perspective
# random_angles %>%
#   group_by(i) %>%
#   mutate(mean_cs = mean(ab_cs)) %>%
#   ggplot(aes(x=ab_cs)) +
#   geom_density() +
#   geom_vline(aes(xintercept=mean_cs), linetype="dashed") +
#   facet_wrap(~i)

# Distribution of cosine similarities taken externally pairwise.
# A and B are different.
random_angles %>%
  filter(i == j) %>%
  ggplot(aes(x=ab_cs)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(ab_cs)), linetype="dashed") +
  labs(x="cos(a, b)",
       title="2. Diagonal sum") ->
  ra_p2

# Distribution of cosine similarities taken internally pairwise.
# A and B are equal.
random_angles %>%
  filter(i > j) %>%
  ggplot(aes(x=abi_cs)) +
  geom_density() +
  geom_vline(aes(xintercept=mean(abi_cs)), linetype="dashed") +
  labs(x="cos(a, b)",
       title="3. Off-diagonal sum") ->
  ra_p3

grid.arrange(ra_p1, ra_p2, ra_p3, ncol=3)






angle_projections %>%
  ggplot(aes(x=wbg_adiff, y=wbg_cuang * wbg_cudnorm)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_smooth(method="lm") +
  labs(x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
       y=TeX("$||D||cos(\\theta_{X_i, A_i - B_i})$"))

angle_projections %>%
  ggplot(aes(y=wbg_pang, x=wbg_adiff, color=wbg_csang * wbg_csdnorm)) + 
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_vline(xintercept=0, linetype="dotted") +
  geom_point() +
  labs(y=TeX("$cos^*(\\theta_{X_i^*, A_i - B_i})$"),
       x=TeX("$cos(\\theta_{X_i, A_i}) - cos(\\theta_{X_i, B_i})$"),
       color=TeX("$||D||cos(\\theta_{X_i, A_i - B_i})$"))

# When comparing two vectors, there is a within-scaling and a between-scaling
# These are linearly related but they have scale bias (more off for larger values)
data.frame(ex_a=scale(weat_white_glove[1,]) - scale(weat_black_glove[1,]),
           ex_b=scale(weat_white_glove[1,] - weat_black_glove[1,])) %>%
  ggplot(aes(x=ex_a, y=ex_b)) +
  geom_smooth() +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype="dashed")

# In the univariate case, the slope is 1, so this problem isn't visible
data.frame(x=weat_white_glove[1,],
           y=as.vector(scale(weat_white_glove[1,]))) %>%
  ggplot(aes(x, y)) +
  geom_point() +
  geom_smooth() +
  geom_abline(slope=1, intercept=0, linetype="dashed")


# What happens as you grow the sum vector
wgs <- t(apply(weat_white_glove, 1, scale2))[sample(1:nrow(weat_white_glove)),]  # Scaled 
wgs_rand_vrare <- t(apply(glove_d[sample(1:400000, 1000),], 1, scale2))  # Scaled 
wgs_rand_rare <- t(apply(glove_d[sample(1:100000, 1000),], 1, scale2))  # Scaled 
wgs_rand_medium <- t(apply(glove_d[sample(1:20000, 1000),], 1, scale2))  # Scaled 
wgs_rand_common <- t(apply(glove_d[sample(1:5000, 1000),], 1, scale2))  # Scaled 

expand.grid(i=2:100, j=1:32) %>%
  rowwise() %>%
  mutate(b_name = weat_white[j],
         dnorm = norm(colSums(wgs_rand_common[1:i,]), "2"),
         sim = lsa::cosine(weat_white_glove[j,], colSums(wgs_rand_common[1:i,]))) %>%
  ggplot(aes(x=dnorm, y=sim, color=i)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(sim)), linetype="dashed") +
  facet_wrap(~b_name)

expand.grid(i=2:100, j=1:32) %>%
  rowwise() %>%
  mutate(b_name = weat_white[j],
         dnorm = norm(colSums(wgs_rand_medium[1:i,]), "2"),
         sim = lsa::cosine(weat_white_glove[j,], colSums(wgs_rand_medium[1:i,]))) %>%
  ggplot(aes(x=dnorm, y=sim, color=i)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(sim)), linetype="dashed") +
  facet_wrap(~b_name)

expand.grid(i=2:100, j=1:32) %>%
  rowwise() %>%
  mutate(b_name = weat_white[j],
         dnorm = norm(colSums(wgs_rand_rare[1:i,]), "2"),
         sim = lsa::cosine(weat_white_glove[j,], colSums(wgs_rand_rare[1:i,]))) %>%
  ggplot(aes(x=dnorm, y=sim, color=i)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(sim)), linetype="dashed") +
  facet_wrap(~b_name)

expand.grid(i=2:100, j=1:32) %>%
  rowwise() %>%
  mutate(b_name = weat_white[j],
         dnorm = norm(colSums(wgs_rand_vrare[1:i,]), "2"),
         sim = lsa::cosine(weat_white_glove[j,], colSums(wgs_rand_vrare[1:i,]))) %>%
  ggplot(aes(x=dnorm, y=sim, color=i)) +
  geom_point() +
  geom_hline(aes(yintercept=mean(sim)), linetype="dashed") +
  facet_wrap(~b_name)


# Plot a single term against an expanding sum vector of random terms from embedding (rarity stratified)
# Note varying X and Y axes.
test_term <- "alex"
expand.grid(i=2:500) %>%
  rowwise() %>%
  mutate(dnorm_common = norm(colSums(wgs_rand_common[1:i,]), "2"),
         dnorm_medium = norm(colSums(wgs_rand_medium[1:i,]), "2"),
         dnorm_rare = norm(colSums(wgs_rand_rare[1:i,]), "2"),
         dnorm_vrare = norm(colSums(wgs_rand_vrare[1:i,]), "2"),
         sim_common = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_common[1:i,])),
         sim_medium = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_medium[1:i,])),
         sim_rare = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_rare[1:i,])),
         sim_vrare = lsa::cosine(glove_d[test_term,], colSums(wgs_rand_vrare[1:i,]))) ->
  term_simd

grid.arrange(term_simd %>%
               ggplot(aes(x=dnorm_common, y=sim_common, color=i)) +
               geom_point() +
               geom_hline(aes(yintercept=mean(sim_common)), linetype="dashed") +
               labs(title=latex2exp::TeX("Components of scalar projection of $v_{Alex}$ onto $D_{rand}$")),
             term_simd %>%
               ggplot(aes(x=dnorm_medium, y=sim_medium, color=i)) +
               geom_point() +
               geom_hline(aes(yintercept=mean(sim_medium)), linetype="dashed") +
               labs(title=latex2exp::TeX("Components of scalar projection of $v_{Alex}$ onto $D_{rand}$")),
             term_simd %>%
               ggplot(aes(x=dnorm_rare, y=sim_rare, color=i)) +
               geom_point() +
               geom_hline(aes(yintercept=mean(sim_rare)), linetype="dashed") +
               labs(title=latex2exp::TeX("Components of scalar projection of $v_{Alex}$ onto $D_{rand}$")),
             term_simd %>%
               ggplot(aes(x=dnorm_vrare, y=sim_vrare, color=i)) +
               geom_point() +
               geom_hline(aes(yintercept=mean(sim_vrare)), linetype="dashed") +
               labs(title=latex2exp::TeX("Components of scalar projection of $v_{Alex}$ onto $D_{rand}$")),
             ncol=2)


# # Run simulation
# set.seed(90210)
# n_samp_min <- 10
# n_samp_max <- 20
# bad_cos_test <- replicate(500, repl_subs(weat_bad, sample(n_samp_min:n_samp_max, 1)), simplify=F)
# bct <- data.frame(do.call(rbind, bad_cos_test)) %>% mutate(value = as.numeric(value))
# bct %>%
#   ggplot(aes(x=value)) +
#   geom_histogram() +
#   facet_grid(t1 ~ t2, scales="free")


## Spectral geometry ##

# Randomize the order; subsample if using
n_samp <- 10
weat_white_glove <- weat_white_glove[sample(1:n_samp),]
weat_black_glove <- weat_black_glove[sample(1:n_samp),]
weat_good_glove <- weat_good_glove[sample(1:n_samp),]
weat_bad_glove <- weat_bad_glove[sample(1:n_samp),]

# Compute SVD to get linear subspaces
# All 4 have different characteristic polynomials
# Scaling the vectors onto the unit ball first does not solve this problem!!
wwg_basis <- svd(weat_white_glove)
wbg_basis <- svd(weat_black_glove)
good_basis <- svd(weat_good_glove)
bad_basis <- svd(weat_bad_glove)

# Compute the matrix that rotates each of these hyperellipses onto the unit ball
wwg_vb <- varimax(wwg_basis$u)
wbg_vb <- varimax(wbg_basis$u)
good_vb <- varimax(good_basis$u)
bad_vb <- varimax(bad_basis$u)


## Shape plots

# First look at the shape of the name space
# Note the white names more heavily weight the larger eigenvalues
name_ev <- rbind(data.frame(x=1:length(wwg_basis$d), y=wwg_basis$d, group="White names"),
                 data.frame(x=1:length(wbg_basis$d), y=wbg_basis$d, group="Black names"))
name_ev %>%
  ggplot(aes(x, y, color=group)) +
  geom_point() +
  geom_line() + 
  labs(x="Eigenvector",
       y="Eigenvalue",
       title="Shape of name subspaces")

# The ratio of eigenvalues in both groups scaled as a proportion of variance provides a sense for
#  how much addition in this domain is scaled by the contribution of each example to the structure of
#  the linear subspace spanned by each set of names (i.e. how differently the A and B subspaces weight each dimension)
# An intuitive way of understanding this is how much you have to stretch or compress your unit ball
# map in each direction to get from one group of words to the other set of words while preserving the meaning of "addition" 
# Both of these elliptical regions can be projected onto the same linear subspace, but they are different projections
d_prop_name <- (wwg_basis$d / sum(wwg_basis$d)) / (wbg_basis$d / sum(wbg_basis$d))
data.frame(x=1:length(wwg_basis$d),
           y=(wwg_basis$d / sum(wwg_basis$d)) / (wbg_basis$d / sum(wbg_basis$d))) %>%
  ggplot(aes(x, y)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_point() +
  labs(x="Eigenvector",
       y="Eigenvector ratio",
       title="Ratio of variance explained by \n dimension between name spaces")

# The ratio of the maximum D_prop to the minimum D_prop describes the maximum distortion.
# When this ratio is 0.5, the eigenvector corresponding to the least influential exemplary pair is half as long
# as the eigenvector for the most influential exemplary pair.
# The scale of the most elliptical 2D cross-section of the hyperellipse
# min(d_prop_name) / max(d_prop_name)

# Now look at the attribute space
# The attribute space is relatively less distorted/stretched than the higher-dimensional name space
attr_ev <- rbind(data.frame(x=1:length(good_basis$d), y=good_basis$d, group="Good words"),
                 data.frame(x=1:length(bad_basis$d), y=bad_basis$d, group="Bad words"))
attr_ev %>%
  ggplot(aes(x, y, color=group)) +
  geom_point() +
  geom_line() + 
  labs(x="Eigenvector",
       y="Eigenvalue",
       title="Shape of attribute subspaces")

# Plot ratio of variance explained by eigenvalue
d_prop_attr <- (good_basis$d / sum(good_basis$d)) / (bad_basis$d / sum(bad_basis$d))
data.frame(x=1:length(good_basis$d),
           y=(good_basis$d / sum(good_basis$d)) / (bad_basis$d / sum(bad_basis$d))) %>%
  ggplot(aes(x, y)) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_point() +
  labs(x="Eigenvector",
       y="Eigenvalue ratio",
       title="Ratio of variance explained by \n dimension between attribute spaces")


## Summary statistics

# For all of the below, look at:
#  - Krz statistic (see above)
#  - Is the basis full rank
#  - Largest/smallest ratio
#  - Changes in cos(a, b) over n for all pairs a, b (how the upper left corner of the matrix changes as it expands)

# TODO: Vary size of word lists (randomly sampling from this word list set; maybe also try doing this in rank order?)
# TODO: Look at two arbitrary sets of words (is this a "placebo test"?)
# TODO: Compare 4 fully random subsets of GloVe


# It can be shown that these have the same dimension, but the eigenvalues are different
# See Krzanowski (1979), "Between-groups comparison of principal components", JASA 74(367): 703-707.
# Compute S matrix (embed in common space); Tr(S) ranges from 0 (orthogonal spaces) to k (coincident spaces)
wg_wb_s <- wwg_basis$u %*% t(wbg_basis$u) %*% wbg_basis$u %*% t(wwg_basis$u)
wbs_basis <- svd(wg_wb_s)
sum(wbs_basis$d)

# The ratio of the maximum D_prop to the minimum D_prop describes the maximum distortion.
# When this ratio is 0.5, the eigenvector corresponding to the least influential exemplary pair is half as long
# as the eigenvector for the most influential exemplary pair.
# The scale of the most elliptical 2D cross-section of the hyperellipse
min(d_prop_name) / max(d_prop_name)
min(d_prop_attr) / max(d_prop_attr)
