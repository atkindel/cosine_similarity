# analogy_test.R
# Implements MSR analogy test from Mikolov et al. to test published embeddings
# Author: Alex Kindel
# Date: 25 March 2022

library(tidyverse)
library(word2vec)
library(lexicon)
library(here)

# Load questions
wr <- read_table2(here("data", "analogy_test", "word_relationship.questions"),
                  col_names = c('t1', 't2', 't3'))
wr.a <- read_table2(here("data", "analogy_test", "word_relationship.answers"),
                    col_names = c('t4'))

# First off, let's explore the analogy questions a bit more
# In reality this test only draws on a very small and selective set of lemmas
# There are 8000 test analogies but only 966 unique tokens and 370 unique lemmas
# This isn't bad for prediction, but has unclear properties for explanation
# One way to think of it is that these lemmas are the conceptual basis that the
#  embedding model space has been trained to express
# Since performance on this task is widely cited as a motivation to use embeddings,
#  published models should perform at least somewhat comparably well on it
wr.tok.ct <- table(unlist(wr))
wr.lem <- tokens_replace(tokens(str_replace(unlist(wr), "'s", "")),
                         pattern = lexicon::hash_lemmas$token,
                         replacement = lexicon::hash_lemmas$lemma)
wr.lem.ct <- table(unlist(wr.lem))

# Note that the tokens in the analogy test are disproportionately drawn from
#  the upper end of the term frequency distribution
# TODO: Term frequency distribution of the tokens

# Construct Levy-Goldberg comparisons + evaluate
# Analogies are of the form t1:t2::t3:?
# Make sure not to consider the question terms as potential answers

# 3CosAdd
wr %>%
  rowwise() %>%
  summarize(word2vec_similarity(as.matrix(nemb[t3,] - nemb[t1,] + nemb[t2,]),
                                as.matrix(nemb)[!rownames(nemb) %in% c(t1, t2, t3),],
                                top_n = 1, type="cosine")) ->
  nemb.analogy.3ca

wr.results.3ca <- cbind(wr, nemb.analogy.3ca)
write(wr.results.3ca$term2, "~/Code/cosine/data/working/nelson_sgns_3ca.txt")
system("~/Code/cosine/data/analogy_test/score.pl ~/Code/cosine/data/working/nelson_sgns_3ca.txt ~/Code/cosine/data/analogy_test/word_relationship.answers")

# PairDirection
# Subtract b off of the target matrix
wr %>%
  rowwise() %>%
  summarize(word2vec_similarity(as.matrix(nemb[t2,] - nemb[t1,]),
                                float::sweep(as.matrix(nemb)[!rownames(nemb) %in% c(t1, t2, t3),], 2, as.numeric(nemb[t3,])),
                                top_n = 1, type="cosine")) ->
  nemb.analogy.pd

wr.results.pd <- cbind(wr, nemb.analogy.pd)
write(wr.results.pd$term2, "~/Code/cosine/data/working/nelson_sgns_pd.txt")
system("~/Code/cosine/data/analogy_test/score.pl ~/Code/cosine/data/working/nelson_sgns_pd.txt ~/Code/cosine/data/analogy_test/word_relationship.answers")

# Original SGNS model is 24.43% accurate on MSR-3CA
# Compare to the Levy & Goldberg SGNS performance benchmark (53.98%)
# Most of the correct performance of the model is in the verb subspace:
#  - Adjectives: 9.53% (45.88%) ~36% difference
#  - Nouns: 18.55% (56.96%) ~38% difference
#  - Verbs: 43.26% (69.90%) ~25% difference

# Note that the above is a little misleading because the given embeddings don't support possessives
# These are about 1/8th of the total analogy set
# But, adding 12.5% (i.e. crediting them all as correct) still doesn't get close.

# Also of interest is inspecting the analogies that relate to terms used in the paper
# This gives us a sense for whether the fixed term lists correspond to what we think
#  they should correspond to in the actual embedding space we have estimated
# Consider e.g. the women and men term lists; the model is only about half as good for
#  the "woman" word list as the "men" word list

women_syn = c('woman', 'women', 'girl', 'girls', 'she', 'her', 'hers', 'herself')
men_syn = c('man', 'men', 'boy', 'boys', 'he', 'him', 'his', 'himself')

wr.results.3ca %>%
  filter(t1 %in% men_syn | t2 %in% men_syn | t3 %in% men_syn) ->
  wr.results.3ca.men  # 23/34 correct (67.64%)

wr.results.3ca %>%
  filter(t1 %in% women_syn | t2 %in% women_syn | t3 %in% women_syn) ->
  wr.results.3ca.women  # 6/18 correct (33.34%)







