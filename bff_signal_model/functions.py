def constant(x, b):
    return linear(x, b, 0)

def linear(x, b, m):
    return m * x + b

def quad(x, b, m, m2):
    return m2 * x ** 2 + m * x + b
