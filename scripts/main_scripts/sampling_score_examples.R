
library(tidyverse)
library(metR)
library(raster)
library(patchwork)

df <- expand.grid(prob=seq(0,1, length.out = 25),
                  uncert=seq(0,1, length.out = 25)) %>% 
  mutate(equal = prob*uncert,
         prob_wt = prob*(uncert^0.3),
         uncert_wt = (prob^0.3)*uncert) %>% 
  pivot_longer(cols = c(3:5), names_to = 'variable', values_to = 'distribs')
head(df)

text_size = 4

## plot equal weighting
p1 <- subset(df, variable == 'equal') %>% 
  ggplot(aes(x=prob, y=uncert, fill = distribs)) +
  geom_tile() +
  geom_vline(xintercept = 0.5, colour = 'red') + geom_hline(yintercept = 0.5, colour = 'red') +
  theme_bw() +
  labs(x='Probability of presence', y='Uncertainty') +
  scale_fill_continuous(type = 'viridis', 
                        name = 'Value of\nDECIDE scores') +
  geom_text(aes(fontface=2), label = 'DECIDE', x = 0.75, y=0.75, colour = 'white', size = text_size) +
  geom_text(aes(fontface=2), label = 'Model-focused', x = 0.25, y=0.75, colour = 'white', size = text_size) +
  geom_text(aes(fontface=2), label = 'User-focused', x = 0.75, y=0.25, colour = 'white', size = text_size) +
  geom_text(aes(fontface=2), label = 'Observations\nunecessary', x = 0.25, y=0.25, colour = 'white', size = text_size) +
  theme(text = element_text(size = 10))
p1


## plot three options
p2 <- ggplot(transform(df, variable = factor(variable, levels = c('prob_wt', 'equal', 'uncert_wt')))) +
  geom_tile(aes(x=prob, y=uncert, fill = distribs)) +
  geom_contour(aes(x=prob, y=uncert, z=distribs), col = 'red') +
  geom_label_contour(aes(x=prob, y=uncert, z=distribs)) +
  facet_grid(~variable,
             labeller = as_labeller(c('equal' = 'Equally weighted',
                                      'uncert_wt' = 'Uncertainty favoured',
                                      'prob_wt' = 'Probability favoured'))) +
  theme_bw() +
  labs(x='Probability of presence', y='Uncertainty') +
  scale_fill_continuous(type = 'viridis', 
                        name = 'Value of\nDECIDE scores') +
  theme(text = element_text(size = 10))



### try a diagram
data <- tibble(x= 1:1000, y= 1:1000)

data %>% 
  ggplot(aes(x, y)) +
  theme_bw() ->
  p
p + inset_element(p1, 0.351, 0.58, 0.61, 0.98,
                  align_to = 'full') +
  inset_element(p2, 0.1, 0.1, 0.9, 0.5) +
  plot_layout(guides = 'collect') +
  

# how does it affect distributions over space?
library(virtualspecies)

# get a raster

# simulate species

# 

## create a landscape
env <- expand.grid(lat=seq(50, 52, length = 25),
                   lon=seq(-3,-1, length = 25)) %>% 
  mutate(value = runif(length(lat))) # how to decide on value?
head(env)

## function?
get_observations <- function(grid, 
                             n_observations,
                             weighting) {
  
  # create ID column 
  grid$ID <- 1:dim(grid)[1]
  
  # generate observations
  obsv <- data.frame(table(sample(1:dim(grid)[1], size = n_observations, replace = T, prob = weighting)))
  
  # match observations to original data
  grid$observations <- obsv$Freq[match(grid$ID, obsv$Var1)]
  
  grid$ID <- NULL
  
  return(grid)
  
}


t <- get_observations(env, 10000, env$value)

head(t)


p1 <- ggplot(data=t, aes(x=lon, y=lat, fill = value)) +
  geom_tile()


p2 <- ggplot(data=t, aes(x=lon, y=lat, fill = observations)) +
  geom_tile()

p1+p2

ggplot(env, aes(x=lon, y=lat, fill = value)) +
  geom_raster(interpolate = F)



r1 <- raster(nrows = 50, ncols = 50, xmn = -3, xmx = -1, ymn = 50, ymx = 52, vals = 0.3)
plot(r1)

rr <- setValues(r1,rnorm(ncell(r1)))
plot(rr)


# get some observations
obsv <- data.frame(table(sample(1:dim(env)[1], size = 10000, replace = T, prob = env$value)))
head(obsv)

# match them
env <- env %>% 
  mutate(obsv = obsv$Freq[match(id, obsv$Var1)])
head(env)

