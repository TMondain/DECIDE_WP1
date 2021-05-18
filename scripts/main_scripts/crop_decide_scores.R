
library(raster)
library(sf)

## load in full decide scores from 
b <- raster(list.files('outputs/decide_score_outputs', full.names = T, pattern = 'butterfly_PA_thinned_10000nAbs')[1])
b_weight <- raster(list.files('outputs/decide_score_outputs', full.names = T, pattern = 'butterfly_PA_thinned_10000nAbs')[3])
m <- raster(list.files('outputs/decide_score_outputs', full.names = T, pattern = 'moth_PA_thinned_10000nAbs')[1])
m_weight <- raster(list.files('outputs/decide_score_outputs', full.names = T, pattern = 'moth_PA_thinned_10000nAbs')[3])

par(mfrow = c(1,2))
plot(m, main = 'moths, unweighted arith mean')
plot(m_weight, main = 'moths weighted presence')

plot(b, main = 'butterflies arith mean')
plot(b, main = 'butterflies weighted presence')

par(mfrow = c(1,1))

# # download map GB
# uk_map <- st_as_sf(getData("GADM", country = "GBR", level = 1, path='Data/environmental_data'))
# uk_map <- st_transform(uk_map, 27700)
# uk_map
# plot(uk_map)
# 
# # remove nrothern ireland
# gb_map <- uk_map[uk_map$NAME_1 != 'Northern Ireland',]
# gb_map
# 
# # check
# plot(st_geometry(gb_map))
# 
# # convert to spatial for use in raster::mask()
# gb_mask <- as_Spatial(gb_map)
# plot(gb_mask[1])
# 
# # moth
# m_gb <- raster::mask(m, gb_mask[1])
# plot(m_gb)
# hist(m_gb)
# 
# # writeRaster(m_gb, 'outputs/decide_score_outputs/moth_GB_decide_score.grd',
# #             format = 'raster', overwrite = T)
# 
# # moth weighted
# m_weight_gb <- raster::mask(m_weight, gb_mask[1])
# plot(m_weight_gb)
# hist(m_weight_gb)
# 
# # writeRaster(m_weight_gb, 'outputs/decide_score_outputs/moth_weighted_prob_pres_GB_decide_score.grd',
# #             format = 'raster', overwrite = T)
# 
# 
# # butterfly
# b_gb <- raster::mask(b, gb_mask[1])
# plot(b_gb)
# hist(b_gb)
# 
# # writeRaster(b_gb, 'outputs/decide_score_outputs/butterfly_GB_decide_score.grd',
# #             format = 'raster', overwrite = T)
# 
# 
# # butterfly weighted
# b_weight_gb <- raster::mask(b_weight, gb_mask[1])
# plot(b_weight_gb)
# hist(b_weight_gb)
# 
# writeRaster(b_weight_gb, 'outputs/decide_score_outputs/butterfly_weighted_prob_pres_GB_decide_score.grd',
#             format = 'raster', overwrite = T)

### read in the files
m_gb <- raster::stack('outputs/decide_score_outputs/moth_GB_decide_score.grd')
m_weight_gb <- raster::stack('outputs/decide_score_outputs/moth_weighted_prob_pres_GB_decide_score.grd')

b_gb <- raster::stack('outputs/decide_score_outputs/butterfly_GB_decide_score.grd')
b_weight_gb <- raster::stack('outputs/decide_score_outputs/butterfly_weighted_prob_pres_GB_decide_score.grd')

par(mfrow = c(1,2))
plot(m_gb, main = 'moths')
plot(b_gb, main = 'butterflies')

plot(m_weight_gb, main = 'moths weighted')
plot(b_weight_gb, main = 'butterflies weighted')


hist(m_gb, main = 'moths')
hist(b_gb, main = 'butterflies')

hist(m_weight_gb, main = 'moths weighted')
hist(b_weight_gb, main = 'butterflies weighted')

par(mfrow = c(1,1))
