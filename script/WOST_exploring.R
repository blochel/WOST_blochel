
# WOST initiation - exploring ----------------------------------------------

library('tidyverse')
library('wader')
library('edenR')
library('tsibble')
library("fpp3")
library('slider')
library('lubridate')
library('mgcv')
library('patchwork')
library('gratia')

# birds -------------------------------------------------------------------


#download_observations("/data")

wost_field_df <- read.csv('data/wost_field_data.cvs',
                          header = TRUE) %>% 
  as_tibble() 

wost_initiate_df <- read.csv('data/stork_initiation.csv',
                             header = TRUE) %>% 
  as_tibble()



WOST_data <- 
  rbind(
    wost_initiate_df %>% 
      filter(initiation %in%
               initiation[str_detect(
                 initiation, "([A-Z]+).*$")]) %>% 
      mutate(month = substr(initiation, 0, 3),
             .before = date_score) ,
    wost_initiate_df %>% 
      filter(initiation %in%
               initiation[str_detect(
                 initiation, "([0-9]+).*$")]) %>% 
      mutate(month = month.abb[month(as.Date(initiation))],
             year = year(as.Date(initiation)),
             .before = date_score)
  ) %>% 
  mutate(new_inititation = 
           match(month, month.abb),
         new_inititation = zoo::as.yearmon(
           paste(year, new_inititation), "%Y %m"),
         .before = month,
         month = factor(month,
                           levels = c('Nov',
                                     'Dec',
                                     'Jan',
                                     'Feb',
                                     'Mar',
                                     'Apr')))


WOST_data %>% 
  ggplot(aes(month))+
geom_histogram(stat="count")
  




# water -------------------------------------------------------------------


depth_3a <- read.csv('data/3as_depth_data.csv',
                     header = TRUE) %>% 
  as_tibble() %>% 
  mutate(time = as.Date(time))

depth_inland <- read.csv('data/inlandenp_depth_data.csv',
                         header = TRUE) %>% 
  as_tibble() %>% 
  mutate(time = as.Date(time))


water_df <- rbind(
depth_3a %>% 
  mutate(area = 'three_A'),
depth_inland %>% 
  mutate(area = 'inland')
) 


water_df %>% 
  ggplot(aes(time, avg_depth, color = area))+
  geom_line()

data <-  dplyr::mutate(water_df, 
                     month=tsibble::yearmonth(time),
                     week=tsibble::yearweek(time))

data2<- data %>% 
  dplyr::group_by(area, week) %>% 
  summarise(mean_depth = mean(avg_depth),
            max_depth = max(avg_depth), 
            min_depth = min(avg_depth)) %>% 
  ungroup() 

data3<- data %>% 
  dplyr::group_by(area, month) %>% 
  summarise(mean_depth = mean(avg_depth),
            max_depth = max(avg_depth), 
            min_depth = min(avg_depth)) %>% 
  ungroup() 




# time decomp -------------------------------------------------------------



# week --------------------------------------------------------------------


# inland ------------------------------------------------------------------


inland_ts <- as_tsibble(data2  %>% 
                         filter(area == 'inland') %>% 
                          select(week,mean_depth), 
                        index = week)
autoplot(inland_ts, mean_depth)




data_ts_inland <- inland_ts %>% 
  mutate(ma_13 = slide_dbl(mean_depth, 
                           mean, 
                           .before = 26, #52 weeks in a year
                           .after = 26,
                           .complete = TRUE #only calculate if you have 6 before and 6 after
  )) 


autoplot(inland_ts, mean_depth) + 
  autolayer(data_ts_inland, ma_13, color = 'blue', size = 2)



add_decomp_inland <- data_ts_inland %>%  model(classical_decomposition(mean_depth, 
                                                         type ='additive')) %>% 
  components()

autoplot(add_decomp_inland) +
  ggtitle('STL - ENP')


stl_output_inland <- data_ts_inland %>% 
  model(STL(mean_depth ~ trend(window = 104) + 
              season(window = 15), 
            robust = TRUE)) %>% 
  components()

autoplot(stl_output_inland) +
  ggtitle('STL - ENP')



# 3A ----------------------------------------------------------------------

threeA_ts <- as_tsibble(data2  %>% 
                          filter(area == 'three_A') %>% 
                          select(week,mean_depth), 
                        index = week)
autoplot(threeA_ts, mean_depth)




data_ts_threeA <- threeA_ts %>% 
  mutate(ma_13 = slide_dbl(mean_depth, 
                           mean, 
                           .before = 26, #52 weeks in a year
                           .after = 26,
                           .complete = TRUE 
  )) 


autoplot(threeA_ts, mean_depth) + 
  autolayer(data_ts_threeA, ma_13, color = 'blue', size = 2)



add_decomp_threeA <- data_ts_threeA %>%  model(classical_decomposition(mean_depth, 
                                                         type ='additive')) %>% 
  components()

autoplot(add_decomp_threeA) +
  ggtitle('STL - 3A')


stl_output_threeA <- data_ts_threeA %>% 
  model(STL(mean_depth ~ trend(window = 104) + 
              season(window = 15), 
            robust = TRUE)) %>% 
  components()

autoplot(stl_output_threeA) +
  ggtitle('STL - 3A')



#first monday of the month 
WOST_weeks <-  dplyr::mutate(WOST_data, 
                       first_of_month = parse_date_time(
                         WOST_data$new_inititation, order = "b Y"),
                       week = tsibble::yearweek(first_of_month), 
                       inittiation = 'yes') %>% 
  select(week, inittiation)




# 3a water + birds ----------------------------------------------------------------------


week_water_birds_3a <- autoplot(stl_output_threeA %>% 
           left_join(WOST_weeks,
                     by = join_by(week)) %>% 
           mutate(week = as.numeric(week))) +
  ggtitle('STL - 3A Scale = Week [trend(104), season(15)]')+
  geom_vline(xintercept = stl_output_threeA %>% 
               left_join(WOST_weeks,
                         by = join_by(week)) %>% 
               filter(inittiation == 'yes') %>% 
               pull(week) %>% as.numeric() ,
             linewidth = 3, color = 'red',
             alpha =0.4)


# enp water + birds---------------------------------------------------------------------


week_water_birds_enp <- autoplot(stl_output_inland %>% 
           left_join(WOST_weeks,
                     by = join_by(week)) %>% 
           mutate(week = as.numeric(week))) +
  ggtitle('STL - ENP Scale = Week [trend(104), season(15)]')+
  geom_vline(xintercept = stl_output_inland %>% 
               left_join(WOST_weeks,
                         by = join_by(week)) %>% 
               filter(inittiation == 'yes') %>% 
               pull(week) %>% as.numeric() ,
             linewidth = 3, color = 'red',
             alpha =0.4)

inland_ts <- as_tsibble(data2  %>% 
                          filter(area == 'inland'), 
                        index = week)



# month --------------------------------------------------------------------


# inland ------------------------------------------------------------------


inland_ts <- as_tsibble(data3  %>% 
                          filter(area == 'inland') %>% 
                          select(month,mean_depth), 
                        index = month)
autoplot(inland_ts, mean_depth)




data_ts_inland <- inland_ts %>% 
  mutate(ma_13 = slide_dbl(mean_depth, 
                           mean, 
                           .before = 7,
                           .after = 7,
                           .complete = TRUE 
  )) 


autoplot(inland_ts, mean_depth) + 
  autolayer(data_ts_inland, ma_13, color = 'blue', size = 2)



add_decomp_inland <- data_ts_inland %>%  model(classical_decomposition(mean_depth, 
                                                                       type ='additive')) %>% 
  components()

autoplot(add_decomp_inland) +
  ggtitle('STL - ENP')


stl_output_inland <- data_ts_inland %>% 
  model(STL(mean_depth ~ trend(window = 23) + 
              season(window = 5), 
            robust = TRUE)) %>% 
  components()

autoplot(stl_output_inland) +
  ggtitle('STL - ENP')



# 3A ----------------------------------------------------------------------

threeA_ts <- as_tsibble(data3  %>% 
                          filter(area == 'three_A') %>% 
                          select(month,mean_depth), 
                        index = month)
autoplot(threeA_ts, mean_depth)




data_ts_threeA <- threeA_ts %>% 
  mutate(ma_13 = slide_dbl(mean_depth, 
                           mean, 
                           .before = 7, 
                           .after = 7,
                           .complete = TRUE 
  )) 


autoplot(threeA_ts, mean_depth) + 
  autolayer(data_ts_threeA, ma_13, color = 'blue', size = 2)



add_decomp_threeA <- data_ts_threeA %>%  model(classical_decomposition(mean_depth, 
                                                                       type ='additive')) %>% 
  components()

autoplot(add_decomp_threeA) +
  ggtitle('STL - 3A')


stl_output_threeA <- data_ts_threeA %>% 
  model(STL(mean_depth ~ trend(window = 23) + 
              season(window = 5), 
            robust = TRUE)) %>% 
  components()

autoplot(stl_output_threeA) +
  ggtitle('STL - 3A')



#first monday of the month 
WOST_months <-  dplyr::mutate(WOST_data, 
                             first_of_month = parse_date_time(
                               WOST_data$new_inititation, order = "b Y"),
                             month = tsibble::yearmonth(first_of_month), 
                             inittiation = 'yes') %>% 
  select(month, inittiation)




# 3a water + birds ----------------------------------------------------------------------


month_water_birds_3a <-autoplot(stl_output_threeA %>% 
           left_join(WOST_months,
                     by = join_by(month)) %>% 
           mutate(month = as.numeric(month))) +
  ggtitle('STL - 3A Scale = month [trend(23), season(5)]')+
  geom_vline(xintercept = stl_output_threeA %>% 
               left_join(WOST_months,
                         by = join_by(month)) %>% 
               filter(inittiation == 'yes') %>% 
               pull(month) %>% as.numeric() ,
             linewidth = 3, color = 'red',
             alpha =0.4)


# enp water + birds ---------------------------------------------------------------------


month_water_birds_enp <-autoplot(stl_output_inland %>% 
           left_join(WOST_months,
                     by = join_by(month)) %>% 
           mutate(month = as.numeric(month))) +
  ggtitle('STL - ENP Scale = month [trend(23), season(5)]')+
  geom_vline(xintercept = stl_output_inland %>% 
               left_join(WOST_months,
                         by = join_by(month)) %>% 
               filter(inittiation == 'yes') %>% 
               pull(month) %>% as.numeric() ,
             linewidth = 3, color = 'red',
             alpha =0.4)



# plots -------------------------------------------------------------------


month_water_birds_3a+
month_water_birds_enp+
week_water_birds_3a+
week_water_birds_enp





# time decomposition ------------------------------------------------------

depth_inland_td <- tsibble::as_tsibble(depth_inland, index=time)

autoplot(depth_inland_td, avg_depth)

water_df_ts <- depth_inland_td %>% 
  mutate(ma_5 = slide_dbl(avg_depth, 
                           mean, 
                           .before = 5,
                           .after = 0,
                           .complete = TRUE), 
         ma_10 = slide_dbl(avg_depth, 
                          mean, 
                          .before = 10,
                          .after = 0,
                          .complete = TRUE),
         msd_5 = slide_dbl(avg_depth, 
                          sd, 
                          .before = 5,
                          .after = 0,
                          .complete = TRUE), 
         msd_10 = slide_dbl(avg_depth, 
                           sd, 
                           .before = 10,
                           .after = 0,
                           .complete = TRUE)) 





# gams  -------------------------------------------------------------------


depth_inland
 #this is the initiation months 

depth_inland_wost <- water_df_ts %>% 
  dplyr::mutate(water_df_ts, 
                first_of_month = parse_date_time(
                  paste(month.abb[month(water_df_ts$time)], year), order = "b Y"),
                month = tsibble::yearmonth(first_of_month)) %>% 
  left_join(WOST_months , 
            join_by(month))


depth_inland_wost <- depth_inland_wost %>% 
  mutate(inittiation = replace_na(inittiation, '0'),
         inittiation = if_else(inittiation == 'yes', '1', inittiation),
         inittiation = as.numeric(inittiation), 
         week = week(time))
  
  



m1 <- gamm(inittiation ~ s(week, bs = "cc") + 
             s(year, k= 10) + s(avg_depth, k = 50),
           data = depth_inland_wost)


summary(m1$gam)
#layout(matrix(1:3, ncol = 3))
#plot(m1$gam, scale = 0)
#layout(1)
draw(m1, residuals = TRUE, rug = FALSE)



m2 <- gam(inittiation ~ s(week, avg_depth, bs = "fs", k = 20) +
            s(ma_10) + #water draw down 10 days
            s(ma_5) + #water draw down 5 days
            #s(msd_10, k = 20)+     
            s(year, k = 15),
               data = depth_inland_wost,
          method = "REML", 
          family = binomial("logit"))
draw(m2, residuals = TRUE, rug = FALSE)


