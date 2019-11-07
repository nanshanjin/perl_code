# Bezier functions

def calculate_coef(p0, p1, p2, p3):
    c = 3*(p1 - p0)
    b = 3*(p2 - p1) -c
    a = p3 - p0 - c - b
    return c, b, a

def Bezier(plist, t):
    # p0 : origin, p1, p2 :control, p3: destination
    p0, p1, p2, p3 = plist
    # calculates the coefficient values
    c, b, a = calculate_coef(p0, p1, p2, p3)
    tsquared = t**2
    tcubic = tsquared*t
    return a*tcubic + b*tsquared + c*t + p0