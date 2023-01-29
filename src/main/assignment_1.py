#import polynomial from numpy for later use in newton-raphson method
from numpy.polynomial import Polynomial

#def functions for absolute and relative error, likely most used functions 
def absolute_error(xi, x):
    return abs(xi-x)

def relative_error(xi, x):
    return absolute_error(xi, x) / abs(xi)

#for questions 1-4
def double_precision_to_decimal(value, chop=-1, rounding=False):
    sign = value[0:1]
    exponent = value[1:12]
    fraction = value[12:]

    s = (-1)**int(sign[0])
    c = 0

    for i in range(len(exponent)):
        c += int(exponent[i]) * 2**(len(exponent)-i-1)

    m = 0

    for i in range(len(fraction)):
        m += int(fraction[i]) * (1/2)**(i+1)

    result = (s) * (2**(c-1023)) * (1+m)
    dot = str(result).find(".")

    if rounding:
        result = float(round(result, chop-dot))
    if chop != -1:
        chopped_result = str(result).replace(".", "")[:chop]
        result = float(chopped_result[:dot] + "." + chopped_result[dot:])

    return result

#prep for question 6
def functionto(x):
    return (x**3) + (4*(x**2)) - 10

#create function for bisection method
def bisection(function, right, left, tol, maxiter):
    i = 0

    while (abs(right - left) > tol and i < maxiter):
        i += 1
        partition = (left + right) / 2
        if (function(left) < 0 and function(partition) > 0) or (function(left) > 0 and function(partition) < 0):
            right = partition
        else:
            left = partition

    return i

#create function for newton raphson method
def nrmethod(previous, tol, maxiter, functionto: Polynomial):
    i = 1
    fncto_derivative = functionto.deriv()

    while i < maxiter:
        if fncto_derivative(previous) == 0:
            return -1
        next_iter = previous - functionto(previous)/fncto_derivative(previous)
        if abs(next_iter - previous) < tol:
            return i

        i += 1
        previous = next_iter
    return i


#first print, trunc to 4
answer1 = double_precision_to_decimal(value="010000000111111010111001" + "0000000000000000000000000000000000000000")
print("%.4f" % answer1, '\n')

#second print, trunc to 2
answer2 = double_precision_to_decimal(value="010000000111111010111001" + "0000000000000000000000000000000000000000", chop=3)
print("%.1f" % answer2, '\n')

#third print, trunc to 2
answer3 = double_precision_to_decimal(value="010000000111111010111001", chop=3, rounding=True)
print("%.1f" % answer3, '\n')

#fourth and fifth print, leave no newline for formatting, no truncation
answer4 = absolute_error(answer1, answer3)
print(answer4)
answer5 = relative_error(answer1, answer3)
print(answer5, '\n')

t = 1
x = 1
while True:
    v = (-1)**t * (x**t/t**3)
    if abs(v) < 10**-4:
        break
    t += 1
print(t-1, '\n')

#print bisection and newton raphson results from functions defined above
print(bisection(function = functionto, right =  7, left =  -4, tol = 10**-4,  maxiter = 100), '\n')
print(nrmethod(-4, 10**-4, 100, Polynomial([-10, 0, 4, 1])), '\n')