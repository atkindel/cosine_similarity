# Replication of "Emotion in political speech"
# Cochrane et al. 2022 Political Communication
# See https://github.com/ccochrane/emotionTranscripts
# 
# Author: Alex Kindel
# Date: 10 April 2022

library(tidyverse)
library(readtext)
library(quanteda)
library(quanteda.textstats)
library(text2vec)
library(gridExtra)
library(lmtest)
library(sandwich)
library(here)

# Load sentiment scores and Hansard data
nemb.d <- read_csv(here("data", "kte2019", "GoogleNews_Embedding.csv"), col_names = FALSE, skip = 1)
nemb.tok <- nemb.d$X1
nemb <- nemb.d[,-1]  # 1.5M X 300d embedding matrix
rownames(nemb) <- nemb.tok
sN <- nrow(nemb)
sp <- ncol(nemb)

# Fixed word lists
affluence <- read_csv(here("data", "kte2019", "word_pairs", "affluence_pairs.csv"), col_names = c("RICH", "POOR"))
gender <- read_csv(here("data", "kte2019", "word_pairs", "gender_pairs.csv"), col_names = c("MASCULINE", "FEMININE"))
race <- read_csv(here("data", "kte2019", "word_pairs", "race_pairs.csv"), col_names = c("BLACK", "WHITE"))

# Load word frequency file
# Corpus is Google Ngram files for year in [2000, 2012]
# TODO

# Load COCA word frequency file
words_219k_m2138 <- read_delim(here("data", "words_219k_m2138.txt"),
                               "\t", escape_double = FALSE, trim_ws = TRUE, skip = 8)

# Link frequencies and compute word vector norms
nt_freq <- data.frame(feature=nemb.tok) %>% rowwise() %>% mutate(snorm = norm(nemb[feature,], "2"))
nt_freq %>%
  left_join(words_219k_m2138, by=c("feature"="word")) %>%
  rename(frequency=freq) ->
  tfn

# Function to plot squared norm as a function of log frequency
# This is Fig. 2 in Arora et al. (2016)
# Also see Eqn. 2.4: log p(w) ~ ||w||/2d - log Z (plus error),
#  where Z = sum exp <w, c>, which is nearly to constant for any choice of context
# This is like the denominator in the softmax function; PMI is a function of a Euclidean
#  distance that has been discounted by the relative maximum distance of this point to
#  everything else (log sum exp over the inner product space at this context) as well as
#  the dimensionality of the inner product space
aroraplot <- function(embm, termfreqs, vocab) {
  termfreqs %>%
    filter(feature %in% vocab) %>%
    rowwise() %>%
    mutate(snorm = norm(embm[feature,], "2"))%>%
    ggplot(aes(x=log(frequency), y=snorm^2)) +
    geom_point() +
    geom_smooth(method="gam")
}

# Dropping one outlier here ("trackback")
# We're totally missing frequencies for the weird ngrams in here like "commenting_policies"
tfn %>%
  filter(!is.na(freq) & snorm <= 10) %>%
  ggplot(aes(x=log(freq), y=snorm^2)) +
  geom_point() +
  geom_smooth(method="gam")

# Squared vector norm is proportional to log term frequency
# The relationship in this model is a little weird for the high frequency terms
# It should be approximately curvilinear but it's quadratic (???)
aroraplot(nemb, tfn, nemb.tok)

# Train a comparable GloVe on this too (takes some time)
# Look at the relationship here as well
# This looks quite a bit more like what we expect to see from Arora et al.
ndf <- fcm(nd)  # Glove uses FCM
ndf@x <- ndf@x + 1  # Use Laplace (add-one) smoothing
nglove <- GlobalVectors$new(rank=100, x_max=100, learning_rate = 0.05)
ngld <- nglove$fit_transform(ndf, n_iter=10, convergence_tol = 0.01)
gemb <- ngld + t(nglove$components)
aroraplot(gemb, nt_freq, rownames(gemb))


# Now let's examine a sample of the inner product space created by this model
# Function to plot inner product manifold given a set of angles/frequencies
plot_angle_manifold <- function(r_angles) {
  r_angles %>%
    ungroup() %>%
    mutate(lfr = log(a$frequency) + log(b$frequency)) %>%
    ggplot(aes(x=nprod, y=ab_cs)) +
    geom_point(aes(color=lfr, size=lfr), alpha=0.9) +
    scale_color_viridis_c() +
    geom_smooth() +
    theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
    labs(x=latex2exp::TeX("$||A||*||B||$"), y=latex2exp::TeX("$cos(A, B)$"),
         color=latex2exp::TeX("$log(p(A)) * log(p(B))$"), size=latex2exp::TeX("$log(p(A)) * log(p(B))$"))
}


# Look at a uniformly random (wrt. term frequency) sample of term pairs
# Default is 5k term pairs
make_angles <- function(embmat, k=5000) {
  sm <- sample(tfn$feature, k*2)
  sA <- sm[1:k]
  sB <- sm[(k+1):(k*2)]
  data.frame(aterm=sA, bterm=sB) %>%
    rowwise() %>%
    mutate(anorm = norm(embmat[aterm,], "2"),
           bnorm = norm(embmat[bterm,], "2"),
           a = tfn[which(tfn$feature == aterm),"frequency"],
           b = tfn[which(tfn$feature == bterm),"frequency"],
           ab_cs = lsa::cosine(as.numeric(embmat[aterm,]), as.numeric(embmat[bterm,]))[1],
           ab_ip = as.numeric(embmat[aterm,]) %*% as.numeric(embmat[bterm,]),
           nprod = anorm * bnorm,
           gK = -1 * (1 + nprod^2 + ab_cs^2)^-2) ->
    random_angles
  
  return(random_angles)
}

ra_nemb <- make_angles(nemb)
plot_angle_manifold(ra_nemb)

# Cosine similarity as a locally linear perspective on the inner product space
# The factorization implies you can look at the manifold from three "sides"
local_linear_decomp <- function(rangs) {
  # Plot some inner product distribution by components. Larger objects in the
  # space have more differentiated IPs, but there are also way fewer of them.
  rangs %>%
    ggplot(aes(y=ab_cs, color=ab_ip, x=nprod, size=ab_ip)) +
    geom_point() +
    scale_color_viridis_c() +
    theme(legend.position="none") ->
    v1
  
  # Cosine similarity differentiates three regions on the inner product space
  # that we can more easily see if we look from the cosine similarity into the
  # manifold. Looking from this perspective helps you see the "twist" structure
  # in A by rotating to the position in which the "flat" structure has a minimal
  # profile and the twist has a maximal profile. The extremal region of the twist
  # (by LNPW) is generally more "flat" than it was from the other perspective.
  # The linear relationship between LNPW and the IP is controlled by cosine similarity.
  rangs %>%
    ggplot(aes(color=ab_cs, y=ab_ip, x=nprod, size=ab_cs)) +
    geom_point() +
    scale_color_viridis_c() +
    theme(legend.position="none") ->
    v2
  
  # Cosine similarity has a linear relationship to the inner product, but only
  # conditional on the local norm product weight, which defines a class of vector pairs
  # with the same CS-IP linear relationship. The distribution of cosine similarities in
  # any given class is generally skewed and tends to have non-zero mean. The variance of
  # the cosine similarity distribution is non-constant and varies systematically with the
  # local norm product weight. One way of thinking about this is that the cosine similarity
  # is a maximally linear perspective on the manifold.
  rangs %>%
    ggplot(aes(x=ab_cs, y=ab_ip, color=nprod, size=nprod)) +
    geom_point() +
    scale_color_viridis_c() +
    theme(legend.position="none") ->
    v3
  
  grid.arrange(v1, v2, v3, ncol=3)
}

local_linear_decomp(ra_nemb)

# 3D visualization makes the hyperbolic-parabolic shape more visible
library(rgl)

# Rescale dimensions of the decomposition
# This makes the curvature easier to see
vid <- view_from_focal_word(nemb, "men")
vid$nprod <- vid$snorm * vid$snorm[1]
vid %<>% filter(similarity < 1)
vid.x <- (vid$nprod - min(vid$nprod)) / (max(vid$nprod) - min(vid$nprod))
vid.y <- (vid$similarity - min(vid$similarity)) / (max(vid$similarity) - min(vid$similarity))
vid.z <- (vid$inner_product - min(vid$inner_product)) / (max(vid$inner_product) - min(vid$inner_product))

# Plot the surface
rgl.open()
rgl.bg(color = "white")
par3d(windowRect = 50 + c(0,0,1000,1000))
rgl.points(vid.x, vid.y, vid.z, color="tomato", size=4, alpha=0.7)
rgl.bbox(color=c("#333377","black"), emission="#333377",
         specular="#3333FF", shininess=5, alpha=0.4)
rgl.lines(c(0, 1), c(0, 0), c(0, 0), color = "red", lwd=2)  # LNPW
rgl.lines(c(0, 0), c(0,1), c(0, 0), color = "blue", lwd=2)  # CS
rgl.lines(c(0, 0), c(0, 0), c(0,1), color = "green", lwd=2)  # IP

# Fit polynomial to this surface with least squares
# We can fit it perfectly since the relationship is fixed
library(rsm)
curvature <- lm(inner_product ~ poly(similarity, nprod, degree=2), data=vid)
persp(curvature, nprod ~ similarity, zlab = "inner_product")
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=30)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=60)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=90)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=120)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=150)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=180)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=210)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=240)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=270)
persp(curvature, nprod ~ similarity, zlab = "inner_product", theta=300)
persp(curvature, nprod ~ similarity, zlab = "inner_product")
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = 10)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = 0)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -10)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -20)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -30)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -40)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -50)
persp(curvature, nprod ~ similarity, zlab = "inner_product", phi = -60)

# Plot cosine simlarity against the log frequency ratio
# I think this one kinda looks like a peacock.
# If you split by LPNW quantiles, you can see the relationship get flat (or close to flat)
ra_nemb %>%
  ungroup() %>%
  mutate(lfr = log(a$frequency) * log(b$frequency)) %>%
  arrange(nprod) %>%
  ggplot(aes(x=lfr, y=ab_cs, color=nprod)) +
  geom_point(size=3) +
  geom_smooth(color="tomato") +
  geom_hline(yintercept=mean(ra_nemb_dim$ab_cs), linetype="dashed") +
  scale_color_viridis_c()

ra_nemb %>%
  ungroup() %>%
  mutate(lfr = log(a$frequency) * log(b$frequency),
         ile = ntile(nprod, 9)) %>%
  arrange(nprod) %>%
  ggplot(aes(x=lfr, y=ab_cs, color=nprod)) +
  geom_point(size=3) +
  geom_smooth(method=MASS::rlm, method.args=list(method="MM"), color="tomato") +
  geom_hline(yintercept=mean(ra_nemb_dim$ab_cs), linetype="dashed") +
  scale_color_viridis_c() +
  facet_wrap(~ile)

# Plot uniformly random cosine similarity sample against log frequency of one component term
# This isn't as informative as the other plots but it's good to know it's not super asymmetric
# it's way more interesting to look at when the A and B sets are highly frequency imbalanced
ra_nemb %>%
  filter(!is.na(a$frequency) & !is.na(b$frequency)) %>%  # Missing term frequency data
  arrange(b$frequency) %>%
  ggplot(aes(x=log(a$frequency), y=ab_cs, color=log(b$frequency))) +
  geom_point() +
  geom_hline(yintercept=mean(ra_nemb$ab_cs), linetype="dashed") +
  geom_smooth() +
  theme(legend.position="bottom") ->
  fp1
ra_nemb %>%
  filter(!is.na(a$frequency) & !is.na(b$frequency)) %>%  # Missing term frequency data
  arrange(a$frequency) %>%
  ggplot(aes(x=log(b$frequency), y=ab_cs, color=log(a$frequency))) +
  geom_point() +
  geom_hline(yintercept=mean(ra_nemb$ab_cs), linetype="dashed") +
  geom_smooth() +
  theme(legend.position="bottom") ->
  fp2
grid.arrange(fp1, fp2, ncol=2)


# Let's look at the Gaussian curvature next
# This shows that the decomposition is non-Euclidean
# It also helps us understand what cosine similarity does to the inner product space
# Mainly what it does is it bends the low-frequency space
# If you scale the LPNW the relative curvature in the rest of the manifold becomes more pronounced
grid.arrange(ra_nemb %>%
  arrange(desc(gK)) %>%
  ggplot(aes(x=nprod, y=ab_cs, color=gK)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, end = 0.95),
ra_nemb %>%
  arrange(desc(gK2)) %>%
  ggplot(aes(x=nprod, y=ab_ip, color=gK)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, end = 0.95),
ra_nemb %>%
  arrange(desc(gK2)) %>%
  ggplot(aes(x=ab_cs, y=ab_ip, color=gK)) +
  geom_point(size=3) +
  scale_color_viridis_c(direction=-1, end = 0.95),
ncol=3)

# Show this in terms of frequency
ra_nemb %>% 
  ggplot(aes(x=log(a$frequency), y=gK, color=ab_cs)) +
  geom_point(size=3) +
  scale_color_viridis_c() +
  theme(legend.position="bottom") -> gkp1
ra_nemb %>% 
  ggplot(aes(x=log(b$frequency), y=gK, color=ab_cs)) +
  geom_point(size=3) +
  scale_color_viridis_c() +
  theme(legend.position="bottom") -> gkp2
grid.arrange(gkp1, gkp2, ncol=2)

# Better plot of curvature-frequency relationship
ra_nemb %>%
  arrange(desc(gK)) %>%
  ggplot(aes(x=log(a$frequency), y=log(b$frequency), color=gK)) +
  geom_point() +
  scale_color_viridis_c(direction=-1)


# Get minimum principal angle
# This quantity is in degrees by default (more easy to think about)
principal_angle <- function(embm, t1, t2, n1, n2, deg=F) {
  sv1 <- svd(embm[c(t1,t2),])
  sv2 <- svd(rbind(embm[t1,]/n1, embm[t2,]/n2))
  dists <- c(geometry::dot(sv1$v[,1], sv2$v[,1]),
             geometry::dot(sv1$v[,2], sv2$v[,1]),
             geometry::dot(sv1$v[,2], sv2$v[,2]),
             geometry::dot(sv1$v[,1], sv2$v[,2]))
  if(deg) {
    return(min(acos(dists)) * (180/pi))
  } else {
    return(min(acos(dists)))
  }
}

# Alternative measure of coincidence
# Krzanowski common space embedding trace statistic
# Sum of squared principal angles between original and norm spaces
krz_trace <- function(embm, t1, t2, n1, n2) {
  d1 <- prcomp(embm[c(t1, t2),])$rotation
  d2 <- prcomp(rbind(embm[t1,]/n1, embm[t2,]/n2))$rotation
  
  krztest <- t(d1) %*% d2 %*% t(d2) %*% d1
  return(sum(svd(krztest)$d))
}

# Or, the minimum angle between vectors in the subspace
krz_angle <- function(embm, t1, t2, n1, n2) {
  d1 <- prcomp(embm[c(t1, t2),])$rotation
  d2 <- prcomp(rbind(embm[t1,]/n1, embm[t2,]/n2))$rotation
  
  krztest <- t(d1) %*% d2 %*% t(d2) %*% d1
  return(acos(max(svd(krztest)$d)))
}

# One helpful way to think about this is in terms of the column space of the embeddings
# Look at the planar basis for the random draws we've done
append_dimensionality <- function(embm, rangs) {
  rangs %>%
    rowwise() %>%
    mutate(vsnorm = norm(embm[aterm,] + embm[bterm,], "2"),
           vdnorm = norm(embm[aterm,] - embm[bterm,], "2"),
           vsnorm.sc = norm(embm[aterm,]/anorm + embm[bterm,]/bnorm, "2"),
           vdnorm.sc = norm(embm[aterm,]/anorm - embm[bterm,]/bnorm, "2"),
           sv1 = svd(embm[c(aterm, bterm),])$d[1],
           sv2 = svd(embm[c(aterm, bterm),])$d[2],
           sv1.sc = svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$d[1],
           sv2.sc = svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$d[2],
           basis.ang = lsa::cosine(svd(embm[c(aterm, bterm),])$v[,1],
                                   svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1]),
           pr.ang = principal_angle(embm, aterm, bterm, anorm, bnorm, deg=T),
           krz = krz_trace(embm, aterm, bterm, anorm, bnorm),
           krza = krz_angle(embm, aterm, bterm, anorm, bnorm),
           d.samedir = sum(sign(embm[aterm,]) * sign(embm[bterm,]) > 0),
           # basis.dir = sum(sign(svd(embm[c(aterm, bterm),])$v[,1]) * sign(svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1]) > 0),
           basis.dir.opp = paste0(which(sign(svd(embm[c(aterm, bterm),])$v[,1]) * sign(svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1]) < 0), collapse="_"),
           basis.dir = ifelse(nchar(basis.dir.opp) == 0, 100, 99 - str_count(basis.dir.opp, "_")),  # Orthant alignment index: 100 - number of opposed dimensions (saves some SVDs)
           basis.dir.2 = ifelse(basis.dir < 50, 100-basis.dir, basis.dir),
           v1 = sv1^2 / (sv1^2 + sv2^2), 
           v2 = sv2^2 / (sv1^2 + sv2^2), 
           v1.sc = sv1.sc^2 / (sv1.sc^2 + sv2.sc^2), 
           v2.sc = sv2.sc^2 / (sv1.sc^2 + sv2.sc^2),
           scaleratio = max(anorm, bnorm)/min(anorm, bnorm)) ->
    rangs_plus
  return(rangs_plus)
}

# Run dimensional analysis for a random subsample (can take a while)
ra_nemb_dim <- append_dimensionality(nemb, ra_nemb %>% ungroup() %>% sample_n(2000))

# Some informative quantities and plots relating to columnar orientation and normalization
# First principal angle: minimum angle between principal vectors of AB and AB'
# Norm deflection score: abs(cos(theta)) between first principal axes of original and normed vector pairs
#  This is the angle between the plane implied by the original vectors and the plane implied by the normalized vectors
#  This value lies approximately on the interval [cos(0), cos(pi/4)] (i.e. parallel to "half orthogonal")
#  So they can be almost parallel or somewhat intersecting
#  We take the absolute value because sometimes we get the normal vector going the other way, but we don't care
#  This is a transformation (hyperbolic centering?) of the first principal angle
#  Note that these are **not** always parallel subspaces! This is super counterintuitive
#  It is very surprising that A and A' are parallel and B and B' are parallel but AB and AB' are not necessarily parallel
#  How is this possible? You have to remember these are in 300 dimensions; 
# Scale ratio: max/min norm{A, B} (how "uneven" or "pointy" is the comparison)
#  There is a clear relationship between the scale ratio and the norm deflection score wrt cos(A, B)
#  In general the normed vectors aren't in the same column space as the original vectors
#  Extremal cosine similarities must have a larger scale ratio when NDS < 1
# Orthant overlap index index: discrete measure on [0,1] of how close they are to being in the same orthant
# 

# First principal angle distribution by scale ratio
# To have a high cosine similarity, the subspaces must be pointy and coincident
ra_nemb_dim %>%
  arrange(desc(abs(ab_cs))) %>%
  ggplot(aes(y=scaleratio, x=pr.ang, color=ab_cs)) +
  geom_point(size=2) +
  scale_color_viridis_c() +
  xlim(0, 90)

# Superimposing the positive and negative NDS regions
# You can see the relationship better this way
ra_nemb_dim %>%
  arrange(desc(abs(ab_cs))) %>%
  ggplot(aes(y=scaleratio, x=abs(basis.ang), color=ab_cs)) +
  geom_point(size=2) +
  geom_vline(xintercept=cos(pi/4), linetype="dashed") +
  scale_color_viridis_c()

# This relationship is self-similar (graph always looks the same no matter where you set the NDS threshold)
ra_nemb_dim %>%
  filter(abs(basis.ang) > 0.95) %>%
  arrange(desc(abs(ab_cs))) %>%
  ggplot(aes(y=scaleratio, x=abs(basis.ang), color=ab_cs)) +
  geom_point(size=2) +
  scale_color_viridis_c() -> bap1
ra_nemb_dim %>%
  filter(abs(basis.ang) > 0.99) %>%
  arrange(desc(abs(ab_cs))) %>%
  ggplot(aes(y=scaleratio, x=abs(basis.ang), color=ab_cs)) +
  geom_point(size=2) +
  scale_color_viridis_c() -> bap2
grid.arrange(bap1, bap2, ncol=2)

# Show NDS on the inner product manifold
# NDS describes Var[cos(A,B)|LPNW(A,B)]
# Smaller NDS values imply cos(A, B) closer to zero
ra_nemb_dim %>%
  arrange(desc(abs(basis.ang))) %>%
  ggplot(aes(x=nprod, y=ab_cs, color=abs(basis.ang), size=-abs(basis.ang))) +
  geom_point() +
  scale_color_viridis_c(direction=-1) +
  scale_size_continuous(range=c(2, 5))

# Using first principal angle
ra_nemb_dim %>%
  arrange(pr.ang) %>%
  ggplot(aes(x=nprod, y=ab_cs, color=pr.ang, size=pr.ang)) +
  geom_point() +
  scale_color_viridis_c() +
  scale_size_continuous(range=c(2, 5))

# Linear relationship between (scaled) norm of summed normalized vectors and spectral norm
# Scaling is linearizing: cosine similarity is perfectly linear in both quantities when rescaled
# In contrast, cosine similarity has "pointiness variance" when not normalized
ra_nemb_dim %>%
  ggplot(aes(x=vsnorm.sc, y=sv1.sc, color=ab_cs)) +
  geom_point() + scale_color_viridis_c() +
  ggtitle("Normalized") + theme(legend.position="bottom") -> scp1
ra_nemb_dim %>%
  ggplot(aes(x=vsnorm, y=sv1, color=ab_cs)) +
  geom_point() + scale_color_viridis_c() +
  ggtitle("Original scale") + theme(legend.position="bottom") -> scp2
grid.arrange(scp1, scp2, ncol=2)

# By component vector % variance explained
# When the space is pointier/longer (?) the cosine similarities become less differentiated, smaller
ra_nemb_dim %>%
  ggplot(aes(x=v1, y=v1.sc, color=ab_cs)) +
  geom_point() +
  scale_color_viridis_c()

# Relationship between first principal angle and NDS
# The principal angle ranges from 0 to 90d; the NDS ranges from 0 to 45d
# When the principal angle is greater than 45d, 
# These are jittered so you can see the concentration of cosines on this manifold
# Low cosines are half-orthogonal; high cosines are close to parallel
# Moderate cosines can be anywhere, sorta
ra_nemb_dim %>%
  arrange(ab_cs) %>%
  mutate(abs.basis.ang.deg = acos(abs(basis.ang)) * (180/pi)) %>%
  ggplot(aes(x=pr.ang, y=abs.basis.ang.deg, color=ab_cs)) +
  geom_jitter(size=2, width=3) +
  geom_vline(xintercept=45, linetype="dashed") +
  scale_color_viridis_c() +
  theme(legend.position="bottom") -> pbr1
ra_nemb_dim %>%
  arrange(desc(ab_cs)) %>%
  mutate(abs.basis.ang.deg = acos(abs(basis.ang)) * (180/pi)) %>%
  ggplot(aes(x=pr.ang, y=abs.basis.ang.deg, color=ab_cs)) +
  geom_jitter(size=2, width=3) +
  geom_vline(xintercept=45, linetype="dashed") +
  scale_color_viridis_c() +
  theme(legend.position="bottom") -> pbr2
grid.arrange(pbr1, pbr2, ncol=2)

# One potentially helpful way to think about it is in terms of the vector space orthants
# Orthants are the n-dimensional equivalent of quadrants; rough descriptor of direction
# Plot the number of directions they have in common against cosine similarity
# Color points by local norm product weight to show the bias
# This has the strong linear relationship we expect
ra_nemb_dim %>%
  arrange(nprod) %>%
  ggplot(aes(y=ab_cs, x=d.samedir, color=nprod)) +
  geom_point(size=3) +
  scale_color_viridis_c()

# Plot orthant overlap index against NDS
# Cosine similarity is lower on average and more variable when the norm-scaled
#  vector subspace and the original-scale vector subspace don't overlap
# Only ~4% pair subspaces of a sample of 2000 even lie in the same orthant of the
#  vector space! Wow
ra_nemb_dim %>%
  arrange(ab_cs) %>%
  ggplot(aes(x=basis.dir.2, y=abs(basis.ang), color=ab_cs)) +
  geom_point(size=3) +
  scale_color_viridis_c()

# Remember that within this orthant there is also some variation in the inner product space
ra_nemb_dim %>%
  filter(basis.dir.2==100) %>%
  ggplot(aes(x=nprod, y=ab_cs)) +
  geom_point()

# Distribution of Krzanowski common subspace embedding trace statistic
# Note concentration at 1
# Note that it appears to be *never* full rank! They have to draw on different parts of the vector space!
ra_nemb_dim %>%
  ggplot(aes(x=krz)) +
  geom_histogram(binwidth=0.01)

# Plot minimum angle between subspaces (Krzanowski method) against scale ratio
# Large cosine similarities tend to be in more square, less coincident subspaces
ra_nemb_dim %>%
  arrange(krza) %>%
  ggplot(aes(y=krza*(180/pi), x=scaleratio, color=ab_cs)) +
  geom_point() +
  geom_smooth(color="gray") +
  scale_color_viridis_c()

# How does this affect arithmetic comparisons of cosine similarities?

# Create many random paired-paired comparisons
# We want to know how the normalization and the frequency bias interact when
#  we compare two paired comparisons
make_angles_geomquad <- function(embm, k=10) {
  sm <- sample(tfn$feature, k*4)
  sA <- sm[1:k]
  sB <- sm[(k+1):(k*2)]
  sC <- sm[(k*2+1):(k*3)]
  sD <- sm[(k*3+1):(k*4)]
  
  data.frame(aterm=sA, bterm=sB, cterm=sC, dterm=sD) %>%
    rowwise() %>%
    mutate(anorm = norm(embm[aterm,], "2"),
           bnorm = norm(embm[bterm,], "2"),
           cnorm = norm(embm[cterm,], "2"),
           dnorm = norm(embm[dterm,], "2"),
           a = tfn[which(tfn$feature == aterm),"frequency"],
           b = tfn[which(tfn$feature == bterm),"frequency"],
           c = tfn[which(tfn$feature == cterm),"frequency"],
           d = tfn[which(tfn$feature == dterm),"frequency"],
           ab_cs = lsa::cosine(as.numeric(embm[aterm,]), as.numeric(embm[bterm,]))[1],
           ab_ip = as.numeric(embm[aterm,]) %*% as.numeric(embm[bterm,]),
           nprod = anorm * bnorm,
           cd_cs = lsa::cosine(as.numeric(embm[cterm,]), as.numeric(embm[dterm,]))[1],
           cd_ip = as.numeric(embm[cterm,]) %*% as.numeric(embm[dterm,]),
           mprod = cnorm * dnorm,
           ab.scaleratio = max(anorm, bnorm)/min(anorm, bnorm),
           cd.scaleratio = max(cnorm, dnorm)/min(cnorm, dnorm),
           ab.pra = principal_angle(embm, aterm, bterm, anorm, bnorm, deg=T),
           cd.pra = principal_angle(embm, cterm, dterm, cnorm, dnorm, deg=T),
           ab.basis.ang = lsa::cosine(svd(embm[c(aterm, bterm),])$v[,1],
                                      svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1]),
           cd.basis.ang = lsa::cosine(svd(embm[c(cterm, dterm),])$v[,1],
                                      svd(rbind(embm[cterm,]/cnorm, embm[dterm,]/dnorm))$v[,1]),
           abcd.ang.sc = lsa::cosine(svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1],
                                     svd(rbind(embm[cterm,]/cnorm, embm[dterm,]/dnorm))$v[,1]),
           abcd.ang = lsa::cosine(svd(rbind(embm[aterm,], embm[bterm,]))$v[,1],
                                      svd(rbind(embm[cterm,], embm[dterm,]))$v[,1])) ->
    random_angles
  
  return(random_angles)
}

gq.nemb <- make_angles_geomquad(nemb, k=4000)

# Plot angle of AB and CD planar bases to each other in original scale (X) and in "cosine scale" (Y)
# Color points by the product of the AB-AB' and CD-CD' angles ("norm deflection score product")
# Split by quantile to show changing relationship; 9iles works well and looks nice, 16iles also useful
# The more parallel the within-comparisons are to their normalized counterparts (NDS product closer to 1),
#  the more the between-comparisons are invariant to normalization. Conversely, when normalization moves
#  the vector pairs into a more different basis, the between-comparison normalized bases tend to be closer
#  to parallel than in their original subspaces, and there is more variance in how close to parallel the
#  corresponding normalized subspaces are. In the more extreme NDSP region, they are less correlated. In 
#  particular, when the NDS product is low, less coincident original subspaces rescale more variably.
gq.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(abs(ab.basis.ang) * abs(cd.basis.ang), 16)) %>%
  arrange(desc(abs(ab.basis.ang) * abs(cd.basis.ang))) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=abs(ab.basis.ang) * abs(cd.basis.ang))) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile)

# Split above plot by AB and BC norm deflection score
gq.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(abs(ab.basis.ang), 16)) %>%
  arrange(desc(abs(ab.basis.ang))) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=abs(ab.basis.ang))) +
  geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile) -> bsp1
gq.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(abs(cd.basis.ang), 16)) %>%
  arrange(desc(abs(cd.basis.ang))) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=abs(cd.basis.ang))) +
  geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile) -> bsp2
grid.arrange(bsp1, bsp2, ncol=2)

# Use untransformed principal angle product instead
# Same pattern (but the quantiles go in the other direction)
gq.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(ab.pra * cd.pra, 16)) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=ab.pra * cd.pra)) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c() + 
  theme(legend.position="bottom") +
  facet_wrap(~ile)

gq.nemb %>%
  ungroup() %>%
  arrange(ab.pra) %>%
  mutate(ile = ntile(ab.pra, 16)) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=ab.pra)) +
  geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile) -> bpp1
gq.nemb %>%
  ungroup() %>%
  arrange(cd.pra) %>%
  mutate(ile = ntile(cd.pra, 16)) %>%
  ggplot(aes(x=abs(abcd.ang), y=abs(abcd.ang.sc), color=cd.pra)) +
  geom_point(size=3) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile) -> bpp2
grid.arrange(bpp1, bpp2, ncol=2)

# Often we are looking at comparisons of the form cos(A, B) @ cos(A, C)
# The shape of these is also interesting
make_angles_geomtri <- function(embm, k=10) {
  sm <- sample(tfn$feature, k*3)
  sA <- sm[1:k]
  sB <- sm[(k+1):(k*2)]
  sC <- sm[(k*2+1):(k*3)]
  
  data.frame(aterm=sA, bterm=sB, cterm=sC) %>%
    rowwise() %>%
    mutate(anorm = norm(embm[aterm,], "2"),
           bnorm = norm(embm[bterm,], "2"),
           cnorm = norm(embm[cterm,], "2"),
           a = tfn[which(tfn$feature == aterm),"frequency"],
           b = tfn[which(tfn$feature == bterm),"frequency"],
           c = tfn[which(tfn$feature == cterm),"frequency"],
           ab_cs = lsa::cosine(as.numeric(embm[aterm,]), as.numeric(embm[bterm,]))[1],
           ab_ip = as.numeric(embm[aterm,]) %*% as.numeric(embm[bterm,]),
           nprod = anorm * bnorm,
           ac_cs = lsa::cosine(as.numeric(embm[aterm,]), as.numeric(embm[cterm,]))[1],
           ac_ip = as.numeric(embm[aterm,]) %*% as.numeric(embm[cterm,]),
           mprod = anorm * cnorm,
           ab.scaleratio = max(anorm, bnorm)/min(anorm, bnorm),
           ac.scaleratio = max(anorm, cnorm)/min(anorm, cnorm),
           ab.basis.ang = lsa::cosine(svd(embm[c(aterm, bterm),])$v[,1],
                                      svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1]),
           ac.basis.ang = lsa::cosine(svd(embm[c(aterm, cterm),])$v[,1],
                                      svd(rbind(embm[aterm,]/anorm, embm[cterm,]/cnorm))$v[,1]),
           ab.pr.ang = principal_angle(embm, aterm, bterm, anorm, bnorm, deg=T),
           ac.pr.ang = principal_angle(embm, aterm, cterm, anorm, cnorm, deg=T),
           abc.ang.sc = lsa::cosine(svd(rbind(embm[aterm,]/anorm, embm[bterm,]/bnorm))$v[,1],
                                    svd(rbind(embm[aterm,]/anorm, embm[cterm,]/cnorm))$v[,1]),
           abc.ang = lsa::cosine(svd(rbind(embm[aterm,], embm[bterm,]))$v[,1],
                                 svd(rbind(embm[aterm,], embm[cterm,]))$v[,1])) ->
    random_angles
  
  return(random_angles)
}

gt.nemb <- make_angles_geomtri(nemb, k=1000)

# The relationship looks a lot worse 
# It gets close to linear when the NDSP is close to 1
# When the NDSP is low, the slope gets close to flat
# Low NDSPs also lead to bimodality (see "split" in panels 1-4)
gt.nemb %>%
  ungroup() %>%
  mutate(ile = ntile(abs(ab.basis.ang) * abs(ac.basis.ang), 12)) %>%
  arrange(desc(abs(ab.basis.ang) * abs(ac.basis.ang))) %>%
  ggplot(aes(x=abs(abc.ang), y=abs(abc.ang.sc), color=abs(ab.basis.ang) * abs(ac.basis.ang))) +
  geom_point(size=3) +
  geom_smooth(method="lm") +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  scale_color_viridis_c(direction=-1) + 
  theme(legend.position="bottom") +
  facet_wrap(~ile)





# Next, look at the distribution of cosine similarities for some arbitrary term
view_from_focal_word <- function(embmat, ft = NA) {
  if(is.data.frame(embmat)) {
    embmat <- as.matrix(embmat)
  }
  if(is.na(ft)) {
    ft <- sample(1:nrow(embmat), 1)
  }

  similmat <- word2vec_similarity(embmat[ft,], embmat, top_n=nrow(embmat), type="cosine")
  
  similmat %>%
    mutate(angle = acos(similarity) * (180/pi)) %>%
    left_join(tfn %>% as.data.frame() %>% dplyr::select(feature, frequency, snorm), by=c("term2"="feature")) %>%
    filter(!is.na(frequency)) %>%  # Drop no-frequency data points
    rowwise() %>%
    mutate(inner_product = embmat[ft,] %*% embmat[term2,]) ->
    test_1term
  return(test_1term)
}

# Plot cosine similarity distribution against log frequency
# Smoothed fits show linear (blue) and generalized additive (red) trends
# Median log frequency indicated by dotted vertical line
# Mean cosine similarity indicated by dashed horizontal line
# In general frequency bias is concentrated in the low-frequency subspace
# The GAM fit in particular is usually almost zero-slope in the high-frequency subspace
# Sometimes it goes back up though (maybe due to high-frequency high-relevance terms, e.g. colors)
# Point color shows vector norm; note somewhat noisy relationship to frequency
plot_focal_view <- function(csims) {
  focal_term <- csims$term2[1]
  csims[-1,] %>%
    ggplot(aes(x=log(frequency), y=similarity, color=snorm)) +
    geom_point() +
    geom_vline(xintercept=median(log(csims$frequency), na.rm=T), linetype="dotted", alpha=0.8) +
    geom_hline(yintercept=mean(csims$similarity, na.rm=T), linetype="dashed", alpha=0.8) +
    geom_smooth(method="gam", color="tomato") +
    geom_smooth(method="lm") +
    scale_color_viridis_c() +
    ggtitle(latex2exp::TeX(sprintf("Frequency-cosine distribution: $\\textit{%s}$", focal_term))) +
    labs(x="Word frequency (log scale)",
         y="Cosine similarity",
         color=latex2exp::TeX("$||w_j||$"))
}

plot_focal_view(view_from_focal_word(nemb))

# From fixed word list
plot_focal_view(view_from_focal_word(nemb, "men"))
plot_focal_view(view_from_focal_word(nemb, "women"))
plot_focal_view(view_from_focal_word(nemb, "white"))
plot_focal_view(view_from_focal_word(nemb, "black"))

# Some cherrypicked/arbitrary ones
# Really they all kind of look like this.
plot_focal_view(view_from_focal_word(nemb, "chews"))
plot_focal_view(view_from_focal_word(nemb, "above"))

# GloVe embedding (note different shape/direction of frequency bias)
plot_focal_view(view_from_focal_word(gemb))

# Local linear decomposition from one term's perspective
local_linear_decomp(view_from_focal_word(nemb, "men") %>%
                      transmute(ab_cs=similarity, ab_ip=inner_product, nprod=snorm*.$snorm[1]) %>%
                      filter(ab_cs < 1))  # We don't care about the self-comparison point




# The analysis is based on arithmetic composites across the term lists
# Let's inspect these directly
# First, modify above functions to work on a composite vetor we provide
view_from_composite <- function(embmat, fv) {
  if(is.data.frame(embmat)) {
    embmat <- as.matrix(embmat)
  }
  
  dotp <- word2vec_similarity(fv, embmat, top_n=nrow(embmat))
  
  word2vec_similarity(fv, embmat, top_n=nrow(embmat), type="cosine") %>%
    mutate(angle = acos(similarity) * (180/pi)) %>%
    left_join(tfn %>% as.data.frame() %>% dplyr::select(feature, frequency, snorm), by=c("term2"="feature")) %>%
    left_join(dotp %>% rename(inner_product=similarity), by="term2") ->
    test_1term
  return(test_1term)
}

plot_composite_view <- function(csims, fv_label, topsel=NA) {
  csims %>%
    ggplot(aes(x=log(frequency), y=similarity, color=snorm)) +
    geom_point() +
    geom_vline(xintercept=median(log(csims$frequency), na.rm=T), linetype="dotted", alpha=0.8) +
    geom_hline(yintercept=mean(csims$similarity, na.rm=T), linetype="dashed", alpha=0.8) +
    geom_smooth(method="gam", color="tomato") +
    geom_smooth(method="lm") +
    scale_color_viridis_c() +
    ggtitle(latex2exp::TeX(sprintf("Frequency-cosine distribution: $\\textit{%s}$", fv_label))) +
    labs(x="Word frequency (log scale)",
         y="Cosine similarity",
         color=latex2exp::TeX("$||w_j||$")) -> plt
  if(is.na(topsel)) {
    print(plt)
  } else {
    print(plt + geom_hline(yintercept=csims[topsel,"similarity"], linetype="dotted"))
  }
}

# For two fixed word lists, get mean of all paired sum vectors
# See Nelson (2021) p. 5
compute_fwl_meanvector <- function(embmat, fwl1, fwl2) {
  expand.grid(fwl1, fwl2) %>%
    rowwise() %>%
    summarize(nemb[which(rownames(nemb) == Var1),] + nemb[which(rownames(nemb) == Var2),]) %>%
    colMeans() ->
    mv_fwl
  return(mv_fwl)
}

# All of the mean vectors exhibit frequency bias as well
mv_male_black <- compute_fwl_meanvector(fwl_male, fwl_black)
mv_female_black <- compute_fwl_meanvector(fwl_female, fwl_black)
mv_male_white <- compute_fwl_meanvector(fwl_male, fwl_white)
mv_female_white <- compute_fwl_meanvector(fwl_female, fwl_white)
plot_composite_view(view_from_composite(nemb, mv_male_black), "male + black (mean vector)")
plot_composite_view(view_from_composite(nemb, mv_female_black), "female + black (mean vector)")
plot_composite_view(view_from_composite(nemb, mv_male_white), "male + white (mean vector)")
plot_composite_view(view_from_composite(nemb, mv_female_white), "female + white (mean vector)")

# It's also useful to look at what we're averaging over
# The biases on each of the component 2sum vectors vary
plot_composite_view(view_from_composite(nemb, as.matrix(nemb["man",] + nemb["black",])), "man + black")
plot_composite_view(view_from_composite(nemb, as.matrix(nemb["his",] + nemb["black",])), "his + black")
plot_composite_view(view_from_composite(nemb, as.matrix(nemb["herself",] + nemb["colored",])), "herself + colored")
plot_composite_view(view_from_composite(nemb, as.matrix(nemb["she",] + nemb["white",])), "she + white")

# Also look at the "social institution" vectors
# The top 50 terms are selected; indicate this with dotted line
# The commented out ylims show the trend in this selected set, roughly
v_polity <- as.matrix(nemb["nation",] + nemb["state",])
v_economy <- as.matrix(nemb["money",])
v_culture <- as.matrix(nemb["culture",])
v_domestic <- as.matrix(nemb["housework",] + nemb["children",])
plot_composite_view(view_from_composite(nemb, v_polity), "nation + state", topsel = 50) #+ ylim(0.638, 1)
plot_composite_view(view_from_composite(nemb, v_economy), "money", topsel = 50) #+ ylim(0.622, 1)
plot_composite_view(view_from_composite(nemb, v_culture), "culture", topsel = 50) #+ ylim(0.665, 1)
plot_composite_view(view_from_composite(nemb, v_domestic), "housework + children", topsel = 50) #+ ylim(0.665, 1)

# How frequency bias affects regression: omitted variable bias in Fig. 3
# Predict difference in means from log term frequency
# Each DiM is a linear regression of the form cos(identity vector, institution vector) ~ identity category
# The frequency bias can be captured by adding various frequency measures to the design matrix

# Construct all 24 comparisons shown in Fig. 3
k <- 50
polity_vmat <- view_from_composite(nemb, v_polity) %>%
  filter(!term2 %in% c("nation", "state") & !is.na(frequency)) %>%
  slice_head(n=k)
polity_words <- polity_vmat$term2
polity_snorm <- polity_vmat$snorm
polity_vmat <- as.matrix(nemb[polity_vmat$term2,])
economy_vmat <- view_from_composite(nemb, v_economy) %>%
  filter(!term2 %in% c("money") & !is.na(frequency)) %>%
  slice_head(n=k)
economy_words <- economy_vmat$term2
economy_snorm <- economy_vmat$snorm
economy_vmat <- as.matrix(nemb[economy_vmat$term2,])
culture_vmat <- view_from_composite(nemb, v_culture) %>%
  filter(!term2 %in% c("culture") & !is.na(frequency)) %>%
  slice_head(n=k)
culture_words <- culture_vmat$term2
culture_snorm <- culture_vmat$snorm
culture_vmat <- as.matrix(nemb[culture_vmat$term2,])
domestic_vmat <- view_from_composite(nemb, v_domestic) %>%
  filter(!term2 %in% c("housework", "children") & !is.na(frequency)) %>%
  slice_head(n=k)
domestic_words <- domestic_vmat$term2
domestic_snorm <- domestic_vmat$snorm
domestic_vmat <- as.matrix(nemb[domestic_vmat$term2,])

cs.male_black_polity <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(polity_vmat[j,], mv_male_black))
cs.male_black_economy <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(economy_vmat[j,], mv_male_black))
cs.male_black_culture <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(culture_vmat[j,], mv_male_black))
cs.male_black_domestic <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(domestic_vmat[j,], mv_male_black))

cs.female_black_polity <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(polity_vmat[j,], mv_female_black))
cs.female_black_economy <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(economy_vmat[j,], mv_female_black))
cs.female_black_culture <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(culture_vmat[j,], mv_female_black))
cs.female_black_domestic <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(domestic_vmat[j,], mv_female_black))

cs.male_white_polity <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(polity_vmat[j,], mv_male_white))
cs.male_white_economy <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(economy_vmat[j,], mv_male_white))
cs.male_white_culture <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(culture_vmat[j,], mv_male_white))
cs.male_white_domestic <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(domestic_vmat[j,], mv_male_white))

cs.female_white_polity <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(polity_vmat[j,], mv_female_white))
cs.female_white_economy <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(economy_vmat[j,], mv_female_white))
cs.female_white_culture <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(culture_vmat[j,], mv_female_white))
cs.female_white_domestic <- sapply(1:nrow(polity_vmat), function(j) lsa::cosine(domestic_vmat[j,], mv_female_white))

make_dim_df <- function(a, b, alab, blab, av, bv, flex, snorm) {
  an <- norm(av, "2")
  bn <- norm(bv, "2")
  rbind.data.frame(data.frame(cs=a, identity=alab, word=flex) %>% 
                     left_join(tfn, by=c("word"="feature")) %>%
                     mutate(nprod = snorm * an),
                   data.frame(cs=b, identity=blab, word=flex) %>%
                     left_join(tfn, by=c("word"="feature")) %>%
                     mutate(nprod = snorm * bn))
}

tsig <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "x")))
}

linear_frequency_bias <- function(dimdf, use.scale=F, use.nprod=T, robust=F, ftest=F, stargazer=F, focal.p=F, fullsumm=T, delta.rsq=F) {
  if(use.scale) {
    dimdf$fbmeasure <- dimdf$snorm^2
  } else if(use.nprod) {
    dimdf$fbmeasure <- 1/dimdf$nprod  
  } else {
    dimdf$fbmeasure <- log(dimdf$frequency)
  }
  ols <- lm(cs ~ identity, data=dimdf)
  ols.freqbias <- lm(cs ~ identity * fbmeasure, data=dimdf)
  if(ftest) {
    anova(ols, ols.freqbias)
  } else if(stargazer) {
    stargazer::stargazer(ols, ols.freqbias)
  } else if(focal.p) {
    summary(polity_mb_fb.freqbias)$coefficients[2,4]
  } else if(fullsumm) {
    if(robust) {
      print(coeftest(ols, vcov = vcovHC(ols, type="HC1")))
      print(coeftest(ols.freqbias, vcov = vcovHC(ols.freqbias, type="HC1")))
    } else {
      print(summary(ols))
      print(summary(ols.freqbias))
      summarize_ovb(dimdf)
    }
  } else if(delta.rsq) {
    data.frame(rsq=c(summary(ols)$r.squared, summary(ols.freqbias)$r.squared),
               model=c("original model", "frequency interaction"))
  } else {
    rbind.data.frame(data.frame(summary(ols)$coefficients, model="original model"),
                     data.frame(summary(ols.freqbias)$coefficients, model="frequency interaction")) %>%
      mutate(sig=tsig(`Pr...t..`))
  }
}

summarize_ovb <- function(dimdf) {
  print(sprintf("Linear correlation of cosine similarity and log word frequency, %s: %f", unique(dimdf$identity)[1],
                cor(dimdf %>% slice_head(n=50) %$% cs, log(dimdf %>% slice_head(n=50) %$% log(frequency)))))
  print(sprintf("Linear correlation of cosine similarity and log word frequency, %s: %f",unique(dimdf$identity)[2],
                cor(dimdf %>% slice_tail(n=50) %$% cs, log(dimdf %>% slice_tail(n=50) %$% log(frequency)))))
  print(sprintf("Covariance of local norm product weight and group (i.e. relative frequency): %f",
                cov(dimdf$identity == unique(dimdf$identity)[1], dimdf$nprod)))
}

# Frequency bias tends to create omitted variable bias because frequency is correlated
#  with the cosine similarity and the relative frequency distributions differ by group
# Adding a measure of frequency into the model tends to increase the standard error on the
#  grouping variable in the difference in mean, and the model fit goes up a lot
# Note that this doesn't always make the effect "go away" per se; the problem is that
#  frequency has a somewhat unpredictable relationship to this problem because it depends
#  on the exact difference in mean log frequency by group and the linear correlation of 
#  the log frequency with the cosine similarity
# In general the reciprocal local norm product weight is the best frequency measure to use
#  because it's a component of cosine similarity, but any functional form will work; the 
#  key thing here is that the LPNW tends to vary by group (frequency bias is relative)

# Black male <-> Black female
linear_frequency_bias(make_dim_df(cs.male_black_polity, cs.female_black_polity, "male_black_polity", "female_black_polity", mv_male_black, mv_female_black, polity_words))
linear_frequency_bias(make_dim_df(cs.male_black_economy, cs.female_black_economy, "male_black_economy", "female_black_economy", mv_male_black, mv_female_black, economy_words))
linear_frequency_bias(make_dim_df(cs.male_black_culture, cs.female_black_culture, "male_black_culture", "female_black_culture", mv_male_black, mv_female_black, culture_words))
linear_frequency_bias(make_dim_df(cs.male_black_domestic, cs.female_black_domestic, "male_black_domestic", "female_black_domestic", mv_male_black, mv_female_black, domestic_words))

# White male <-> White female
linear_frequency_bias(make_dim_df(cs.male_white_polity, cs.female_white_polity, "male_white_polity", "female_white_polity", mv_male_white, mv_female_white, polity_words))
linear_frequency_bias(make_dim_df(cs.male_white_economy, cs.female_white_economy, "male_white_economy", "female_white_economy", mv_male_white, mv_female_white, economy_words))
linear_frequency_bias(make_dim_df(cs.male_white_culture, cs.female_white_culture, "male_white_culture", "female_white_culture", mv_male_white, mv_female_white, culture_words))
linear_frequency_bias(make_dim_df(cs.male_white_domestic, cs.female_white_domestic, "male_white_domestic", "female_white_domestic", mv_male_white, mv_female_white, domestic_words))

# White female <-> Black female
linear_frequency_bias(make_dim_df(cs.female_white_polity, cs.female_black_polity, "female_white_polity", "female_black_polity", mv_female_white, mv_female_black, polity_words))
linear_frequency_bias(make_dim_df(cs.female_white_economy, cs.female_black_economy, "female_white_economy", "female_black_economy", mv_female_white, mv_female_black, economy_words))
linear_frequency_bias(make_dim_df(cs.female_white_culture, cs.female_black_culture, "female_white_culture", "female_black_culture", mv_female_white, mv_female_black, culture_words))
linear_frequency_bias(make_dim_df(cs.female_white_domestic, cs.female_black_domestic, "female_white_domestic", "female_black_domestic", mv_female_white, mv_female_black, domestic_words))

# Black male <-> White male
linear_frequency_bias(make_dim_df(cs.male_black_polity, cs.male_white_polity, "male_black_polity", "male_white_polity", mv_male_black, mv_male_white, polity_words))
linear_frequency_bias(make_dim_df(cs.male_black_economy, cs.male_white_economy, "male_black_economy", "male_white_economy", mv_male_black, mv_male_white, economy_words))
linear_frequency_bias(make_dim_df(cs.male_black_culture, cs.male_white_culture, "male_black_culture", "male_white_culture", mv_male_black, mv_male_white, culture_words))
linear_frequency_bias(make_dim_df(cs.male_black_domestic, cs.male_white_domestic, "male_black_domestic", "male_white_domestic", mv_male_black, mv_male_white, domestic_words))

# White male <-> Black female
linear_frequency_bias(make_dim_df(cs.male_white_polity, cs.female_black_polity, "male_white_polity", "female_black_polity", mv_male_white, mv_female_black, polity_words))
linear_frequency_bias(make_dim_df(cs.male_white_economy, cs.female_black_economy, "male_white_economy", "female_black_economy", mv_male_white, mv_female_black, economy_words))
linear_frequency_bias(make_dim_df(cs.male_white_culture, cs.female_black_culture, "male_white_culture", "female_black_culture", mv_male_white, mv_female_black, culture_words))
linear_frequency_bias(make_dim_df(cs.male_white_domestic, cs.female_black_domestic, "male_white_domestic", "female_black_domestic", mv_male_white, mv_female_black, domestic_words))

# White female <-> Black male
linear_frequency_bias(make_dim_df(cs.female_white_polity, cs.male_black_polity, "female_white_polity", "male_black_polity", mv_female_white, mv_male_black, polity_words))
linear_frequency_bias(make_dim_df(cs.female_white_economy, cs.male_black_economy, "female_white_economy", "male_black_economy", mv_female_white, mv_male_black, economy_words))
linear_frequency_bias(make_dim_df(cs.female_white_culture, cs.male_black_culture, "female_white_culture", "male_black_culture", mv_female_white, mv_male_black, culture_words))
linear_frequency_bias(make_dim_df(cs.female_white_domestic, cs.male_black_domestic, "female_white_domestic", "male_black_domestic", mv_female_white, mv_male_black, domestic_words))





