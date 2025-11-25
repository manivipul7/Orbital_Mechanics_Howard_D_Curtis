import numpy as np

def bisect(fun, xl, xu):
    """
    This function evaluates a root of a function using the bisection method.

    Parameters:
    fun -- the function whose root is being found
    xl -- low end of the interval containing the root
    xu -- upper end of the interval containing the root

    Returns:
    root -- the computed root
    """
    tol = 1.e-6
    n = int(np.ceil(np.log(np.abs(xu - xl) / tol) / np.log(2)))

    for i in range(n):
        xm = (xl + xu) / 2.0
        fxl = fun(xl)
        fxm = fun(xm)
        if fxl * fxm > 0:
            xl = xm
        else:
            xu = xm

    root = xm
    return root

# Example usage:
if __name__ == "__main__":
    # Define a sample function, e.g., f(x) = x^2 - 4
    def f(x):
        return x**2 - 4

    # Set the interval [xl, xu]
    xl = 1
    xu = 3

    # Call the bisect function
    root = bisect(f, xl, xu)
    print(f"The root is: {root}")
