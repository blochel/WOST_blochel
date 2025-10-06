library("gratia")
library("mgcv")
library('tidyverse')

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
             i.initiation = 2,
             .before = date_score) ,
    wost_initiate_df %>% 
      filter(initiation %in%
               initiation[str_detect(
                 initiation, "([0-9]+).*$")]) %>% 
      mutate(month = month.abb[month(as.Date(initiation))],
             year = year(as.Date(initiation)),
             i.initiation = 2, 
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




WOST_months <- dplyr::mutate(WOST_data, 
              first_of_month = parse_date_time(
                WOST_data$new_inititation, order = "b Y"),
              month = tsibble::yearmonth(first_of_month), 
              inittiation = 'yes')



depth_inland_wost <- data %>% 
  left_join(WOST_months , 
            join_by(month)) %>% 
  mutate(months = month(time), 
         weeks = week(time))

depth_inland_wost$i.initiation[is.na(depth_inland_wost$i.initiation)] <- 1




wost_inland <- depth_inland_wost %>% 
  filter(area == 'inland') %>% 
  dplyr::select(-c(year.x, year.y, initiation, new_inititation, 
                   date_score, days_past_nov_1, region, colony, 
                   notes, first_of_month, inittiation)) %>% 
  mutate(i.initiation = as.integer(i.initiation),
         doy = as.numeric(strftime(time, format = "%j")),
         f.year = as.factor(year(time)),
         bi_initiation = if_else(i.initiation == 1, 0, 1)) %>% 
  group_by(month) %>% 
  mutate(max_m_depth = mean(avg_depth, na.omit = TRUE)) %>% 
  ungroup() %>% 
  group_by(bi_initiation) %>% 
  mutate( drawd_1 = lag(avg_depth, 1) - avg_depth,
          drawd_2 = lag(avg_depth, 2) - avg_depth,
          drawd_3 = lag(avg_depth, 3) - avg_depth,
          drawd_4 = lag(avg_depth, 4) - avg_depth,
          drawd_7 = lag(avg_depth, 7) - avg_depth,
          drawd_14 = lag(avg_depth, 14) - avg_depth) %>% 
  ungroup() 





# The reproductive period takes 120 to 130 days (approximately four months). 
# The breeding season refers to the period between the time the female lays 
# the first egg and the time the first nestling leaves the nest.



# gammmtime# gammming  ---------------------------------------------------------------





m_cat_data1 <- gam(i.initiation ~ s(doy, drawd_7,  k = 60)+ #great
                   #s(bi_initiation, drawd_14) +
                   #s(weeks) + 
                   s(months, bs = 'cc', k = 12), #good 
                   data = wost_inland,              
                   family = ocat(R=3),       
                   method = "REML"
)    


summary(m_cat_data1)
draw(m_cat_data1, residuals = TRUE, rug = FALSE) #look these up under GRATIA
draw(m_cat_data1, unconditional = TRUE, parametric = FALSE) #look these up 

gam.check(m_cat_data1)




pred_aug_1 <- data.frame(doy = yday(ymd('2030-01-20')), 
                         #weeks = week(ymd('2010-01-20')),
                         months = month(ymd('2030-01-20')), 
                         drawd_7 = -2)
predict(m_cat_data1, newdata = pred_aug_1, type = 'response', se=TRUE)

# data slice over doy/months(!) holding other covariates at median values
ds1.1 <- data_slice(m_cat_data1, months = evenly(months, n = 100))
fv1.1 <- fitted_values(m_cat_data1, data = ds1.1)

fv1.1 %>% 
  ggplot(aes(x = months, y = .fitted, group = .category)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .category),
              alpha = 0.2) +
  geom_line(aes(colour = .category))


#M 1.2 days and drawdown 

m_cat_data1.2 <- gam(i.initiation ~ s(doy, drawd_7,  k = 50), #+ #great
                   #s(bi_initiation, drawd_14) +
                   #s(weeks) + 
                   #s(months, bs = 'cc'), 
                   data = wost_inland,              
                   family = ocat(R=3),       
                   method = "REML"
)    


summary(m_cat_data1.2)
draw(m_cat_data1.2, residuals = TRUE, rug = FALSE)

k.check(m_cat_data1.2)
draw(m_cat_data1.2, unconditional = TRUE, parametric = FALSE)



pred_aug_1 <- data.frame(doy = yday(ymd('2030-01-20')), 
                         #weeks = week(ymd('2010-01-20')),
                         #months = month(ymd('2030-01-20')), 
                         drawd_7 = -2)
predict(m_cat_data1.2, newdata = pred_aug_1, type = 'response', se=TRUE)

# data slice over doy/months(!) holding other covariates at median values
ds1.2 <- data_slice(m_cat_data1.2, doy = evenly(doy, n = 100))
fv1.2 <- fitted_values(m_cat_data1.2, data = ds1.2)

fv1.2 %>% 
  ggplot(aes(x = doy, y = .fitted, group = .category)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .category),
              alpha = 0.2) +
  geom_line(aes(colour = .category))



# adding year -------------------------------------------------------------



m_cat_data2 <- bam(i.initiation ~  s(months) + s(drawd_14) + 
                     s(months, f.year, bs = 'sz'), 
        data = wost_inland,              
        family = ocat(R=3),    
        discrete = TRUE, 
        nthreads = 6,
        method = "fREML"
)    

summary(m_cat_data2)
draw(m_cat_data2, residuals = TRUE, rug = FALSE)


# data slice over months/weeks(!) holding other covariates at median values
ds2 <- data_slice(m_cat_data2, months = evenly(months, n = 100),
                  drawd_14 = evenly(drawd_14, n = 100))
fv2 <- fitted_values(m_cat_data2, data = ds2)

fv2 %>% 
  ggplot(aes(x = months, y = .fitted, group = .category)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, fill = .category),
              alpha = 0.2) +
  geom_line(aes(colour = .category))




# testing binomial data ---------------------------------------------------


m_bi_data <- bam(bi_initiation ~  s(weeks) + s(drawd_7) + 
            te(weeks, drawd_7, bs = c("cc",'re')), #what is this cc and fs?
          data = wost_inland,       
          method = "REML",               
          family = binomial("logit"))   
summary(m_bi_data)
draw(m_bi_data, residuals = TRUE, rug = FALSE)



ds3 <- data_slice(m_bi_data, weeks = evenly(weeks, n = 100),
                  drawd_7 = evenly(drawd_7, n = 100))
fv3 <- fitted_values(m_bi_data, data = ds3)


fv3 %>% 
  group_by(weeks) %>% 
  summarise(weeks = weeks,
            .fitted = .fitted,
            .lower_ci = .lower_ci,
            .upper_ci = .upper_ci) %>% 
  ggplot(aes(x = weeks, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              alpha = 0.2) 


fv3 %>% 
  group_by(weeks) %>% 
  summarise(weeks = mean(weeks),
            .fitted = mean(.fitted),
            .lower_ci = mean(.lower_ci),
            .upper_ci = mean(.upper_ci)) %>% 
  ggplot(aes(x = weeks, y = .fitted)) +
  geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci),
              alpha = 0.2) +
  geom_line(color = 'blue', linewidth = 2)




#structure test data around the nesting Nov - April. 
#






test_data <- wost_inland %>% 
  filter( year(time) <= 2023)

# Predictions with standard errors
pred <- predict(m_bi_data, newdata = test_data, type = "response", se.fit = TRUE)
fit_vals <- pred$fit
se_vals  <- pred$se.fit
lower <- fit_vals - 1.96 * se_vals
upper <- fit_vals + 1.96 * se_vals





# Quick plot
plot(test_data$time, test_data$bi_initiation, pch = 16, col = "blue",
     xlim = c(min(test_data$time), max(test_data$time)), ylim = range(test_data$bi_initiation))
lines(test_data$time, fit_vals, col = "red", lwd = 2)
lines(test_data$time, lower, col = "grey", lty = 2)
lines(test_data$time, upper, col = "grey", lty = 2)



