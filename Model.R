##List of Common Machine Learning Algorithms

## 1. Linear Regression
## Y= a *X + b.


#Load Train and Test datasets
#Identify feature and response variable(s) and values must be numeric and numpy arrays
x_train <- data.frame(x_train= c(2, 4, 5, 6, 7, 3, 4, 6)) ##input_variables_values_training_datasets
y_train <- data.frame(y_train = c(4, 7, 11, 13, 14, 6, 8.5, 13.5)) ##target_variables_values_training_datasets
x_test <- data.frame(x_train= c(1, 20, 30, 40, 50, 60, 4, 6)) ##input_variables_values_test_datasets
x <- cbind(x_train,y_train)
# Train the model using the training sets and check score
linear <- lm(x$y_train ~ x_train, data = x)
summary(linear)
#Predict Output
predicted= predict(linear,newdata = x_test) 




## 2. Logistic Regression
#  it predicts the probability of occurrence of an event by fitting data to a logit function. 
# Hence, it is also known as logit regression. Since, it predicts the probability, its output 
# values lies between 0 and 1 (as expected).

# odds= p/ (1-p) = probability of event occurrence / probability of not event occurrence
# ln(odds) = ln(p/(1-p))
# logit(p) = ln(p/(1-p)) = b0+b1X1+b2X2+b3X3....+bkXk


x_train <- data.frame(x_train = c(4, 7, 11, 13, 14, 6, 8.5, 13.5))
y_train <- data.frame(y_train = c(0, 0, 1, 1, 1, 0, 0, 1)) ##target_variables_values_training_datasets
x <- cbind(x_train,y_train)
x_test <- data.frame(x_train= c(1, 20, 30, 40, 50, 60, 4, 6)) ##input_variables_values_test_datasets

# Train the model using the training sets and check score
logistic <- glm(y_train ~ x_train , data = x,family='binomial')
summary(logistic)
#Predict Output
predicted= predict(logistic,newdata = x_test)
predicted