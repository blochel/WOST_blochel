
# Forecasting with binomial data YE/NO ----------------

library("gratia")
library("mgcv")
library('tidyverse')




# for ggplots -------------------------------------------------------------

theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(axis.line.x.bottom = element_line(colour = "black",
                                                    size = 1),
                  axis.line.y.left = element_line(colour = "black",
                                                  size = 1)))


# data download -----------------------------------------------------------


wost_field_df <- read.csv('data/wost_field_data.cvs',
                          header = TRUE) %>% 
  as_tibble() 

wost_initiate_df <- read.csv('data/stork_initiation.csv',
                             header = TRUE) %>% 
  as_tibble()



# wost data ---------------------------------------------------------------



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
      mutate(month = month.abb[month(as.Date(initiation))],                    #make this week!
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



# wost and water ----------------------------------------------------------



wost_inland <- depth_inland_wost %>% 
  filter(area == 'inland') %>% 
  dplyr::select(-c(year.x, year.y, initiation, new_inititation, 
                   date_score, days_past_nov_1, region, colony, 
                   notes, first_of_month, inittiation)) %>% 
  mutate(i.initiation = as.integer(i.initiation),
         doy = as.numeric(strftime(time, format = "%j")),
         f.year = as.factor(year(time)),
         bi_initiation = if_else(i.initiation == 1, 0, 1)) %>% 
  group_by(week) %>% 
  mutate(max_week_depth = mean(avg_depth, na.omit = TRUE)) %>% 
  ungroup() %>% 
  group_by(bi_initiation) %>% 
  mutate( drawd_1 = lag(avg_depth, 1) - avg_depth,
          drawd_2 = lag(avg_depth, 2) - avg_depth,
          drawd_3 = lag(avg_depth, 3) - avg_depth,
          drawd_4 = lag(avg_depth, 4) - avg_depth,
          drawd_7 = lag(avg_depth, 7) - avg_depth,
          drawd_14 = lag(avg_depth, 14) - avg_depth) %>% 
  mutate(year_depth_lag = lag(max_week_depth, 183)-max_week_depth) %>% 
  ungroup() 


plot(wost_inland$time, wost_inland$year_depth_lag)

plot(wost_inland$time, wost_inland$bi_initiation)

# Split data to train and test sets ---------------------------------------

train_prop <- 0.95
wost_train <- wost_inland[1:ceiling(train_prop * nrow(wost_inland)),]
wost_test <- wost_inland[(ceiling(train_prop * nrow(wost_inland)) + 1):nrow(wost_inland),]



plot(wost_train$time, wost_train$bi_initiation, 
     main = 'training data')
plot(wost_test$time, wost_test$bi_initiation, 
     main = 'test data')


# gamming -----------------------------------------------------------------




m1_wost_bi <- bam(bi_initiation ~  
                   s(weeks) +
                   s(drawd_7) + 
                   s(year_depth_lag)+
                    te(year_depth_lag, drawd_7, bs = c("re",'re'))+
                    te(weeks, drawd_7, bs = c("cc",'re'))
                  , 
                 data = wost_train,       
                 method = "REML",               
                 family = binomial("logit"))   
summary(m1_wost_bi)
draw(m1_wost_bi, residuals = TRUE, rug = FALSE)



#ds3 <- data_slice(m1_wost_bi, weeks = evenly(weeks, n = 1000),
#                  drawd_7 = evenly(drawd_7, n = 1000))
#fv3 <- fitted_values(m1_wost_bi, data = ds3)


plot(wost_train$time, wost_train$bi_initiation, 
     main = 'training data')


#not working .... 
#plot(wost_train$time, wost_train$bi_initiation,
#     xlab = "date", ylab = "nesting prob")
#lines(wost_train$time, m1_wost_bi$fitted.values, col = "red")


test_predictions <- predict(m1_wost_bi, newdata = wost_test, se.fit = TRUE) 
plot(wost_test$time, wost_test$bi_initiation,
     xlab = "date", ylab = "nesting prob")
lines(wost_test$time, test_predictions$fit, col = "red")


data.frame(date = wost_test$time, 
           model_prediction = test_predictions$fit) %>% 
  filter(model_prediction >= 0.5) 


WOST_data %>% filter(year >= 2023)

